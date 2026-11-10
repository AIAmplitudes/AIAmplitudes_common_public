"""Relation counts follow selected data and share the relation-loading cache."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import aiamplitudes_common_public as public
from aiamplitudes_common_public import fbspaces as fb


class RelationCountTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.sources = public.frontspace_file, dict(public._backspace_files)
        public.clear_basis_cache()

    def tearDown(self):
        public.frontspace_file, backfiles = self.sources
        public._backspace_files.update(backfiles)
        public.clear_basis_cache()

    def relation_file(self, filename, equations, *, space='front', weight=2,
                      label='customrels'):
        path = self.root / filename
        path.parent.mkdir(parents=True, exist_ok=True)
        header = 'init' if space == 'front' else 'fin'
        names = ','.join(['unused'] * (weight - 1) + [label])
        path.write_text(f'{header}_rels := [{names}]:\n'
                        f'{header}_list[{weight}] := customindep :\n'
                        f'{label} := [{",".join(equations)}]:\n')
        return path

    def test_legacy_sections_without_metadata(self):
        path = self.root / 'legacy'
        path.write_text(
            '# idoublerels999 := [SB(xx)=0]:\n'
            'idoubleindep2 := [SB(aa),SB(ab)]:\n'
            'idoublerels17 := [SB(ac)=SB(aa),SB(ad)=0]:\n'
            'doubleindep1 := [E(aa)]:\n'
            'doublerels12REF := [E(xx)=0]:\n'
            'doublerels12 := [E(ab)=E(aa)]:\n')
        public.configure_basis(front_rest=path, back_rest=path)
        self.assertEqual(public.basis_size(2, 'front'), 2)
        self.assertEqual(public.basis_size(2, 'back'), 1)
        self.assertEqual(public.relation_count(2, 'front'), 2)
        self.assertEqual(public.relation_count(2, 'back'), 1)

    def test_explicit_metadata_overrides_legacy_labels(self):
        path = self.relation_file('metadata', ['SB(ab)=0'])
        with path.open('a') as handle:
            handle.write('idoubleindep9 := [SB(aa)]:\n'
                         'idoublerels9 := [SB(ac)=0,SB(ad)=0]:\n')
        rels, basis = fb._parse_file_header(str(path), 'init')
        self.assertEqual(rels[2], 'customrels')
        self.assertEqual(basis[2], 'customindep')
        self.assertEqual(public.relation_count(2, prefix=path), 1)

    def test_front_load_then_count_reuses_parsed_data(self):
        path = self.relation_file('front', ['SB(ab)=SB(aa)', 'SB(ac)=0'])
        public.configure_basis(front_rest=path)
        with patch.object(fb, 'getFrel_eqs', wraps=fb.getFrel_eqs) as reader:
            relations = public.fr_rels(2)
            self.assertEqual(relations, {'ab': {'aa': 1}, 'ac': {None: 0}})
            for _ in range(5):
                self.assertEqual(fb.NF_rels[2], 2)
                self.assertEqual(public.relation_count(2), 2)
                self.assertIs(public.fr_rels(2), relations)
                self.assertIs(fb.get_frels(2, path.parent, path.name), relations)
            self.assertEqual(reader.call_count, 1)
        self.assertEqual(list(fb.relation_counts.values()), [2])

    def test_back_count_then_load_reuses_parsed_data(self):
        path = self.relation_file('back', ['E(ab)=E(aa)', 'E(ac)=0'], space='back')
        public.configure_basis(back_rest=path)
        with patch.object(fb, 'getBrel_eqs', wraps=fb.getBrel_eqs) as reader:
            self.assertEqual(fb.NB_rels[2], 2)
            relations = public.br_rels(2)
            for _ in range(5):
                self.assertEqual(public.relation_count(2, 'back'), 2)
                self.assertIs(public.br_rels(2), relations)
            self.assertEqual(reader.call_count, 1)

    def test_switching_selected_files_updates_existing_count_views(self):
        first = self.relation_file('one', ['SB(ab)=0'])
        second = self.relation_file('two', ['SB(ab)=0', 'SB(ac)=0'])
        front_counts = fb.NF_rels
        for path, expected in ((first, 1), (second, 2), (first, 1)):
            public.configure_basis(front_rest=path)
            self.assertEqual(front_counts[2], expected)
            self.assertEqual(len(public.fr_rels(2)), expected)

    def test_explicit_sources_and_directories_do_not_change_defaults(self):
        first = self.relation_file('first/relations', ['SB(ab)=0'])
        second = self.relation_file('second/relations', ['SB(ab)=0', 'SB(ac)=0'])
        public.configure_basis(front_rest=first)
        for path, count in ((first, 1), (second, 2), (first, 1)):
            self.assertEqual(public.relation_count(2, prefix=path), count)
            self.assertEqual(public.relation_count(2, mydir=path.parent,
                                                   prefix=path.name), count)
            self.assertEqual(len(public.fr_rels(2, prefix=path)), count)
        self.assertEqual(fb.NF_rels[2], 1)
        self.assertEqual(len(fb.relation_counts), 2)

    def test_back_operator_traces_use_their_selected_files(self):
        phi2 = self.relation_file('phi2', ['E(ab)=0'], space='back')
        phi3 = self.relation_file('phi3', ['E(ab)=0', 'E(ac)=0'], space='back')
        public.configure_basis(back_rest=phi2)
        public.configure_basis(back_rest=phi3, optrace='phi3')
        self.assertEqual(fb.NB_rels[2], 1)
        self.assertEqual(public.relation_count(2, 'back', optrace='phi3'), 2)
        self.assertEqual(len(public.br_rels(2, optrace='phi3')), 2)

    def test_front_and_back_sections_in_one_file_have_separate_counts(self):
        front = self.relation_file('front', ['SB(ab)=0'])
        back = self.relation_file('back', ['E(ab)=0', 'E(ac)=0'], space='back')
        combined = self.root / 'combined'
        combined.write_text(front.read_text().replace('customrels', 'frontrels') +
                            back.read_text().replace('customrels', 'backrels'))
        self.assertEqual(public.relation_count(2, 'front', prefix=combined), 1)
        self.assertEqual(public.relation_count(2, 'back', prefix=combined), 2)

    def test_clear_reloads_relations_counts_and_headers_after_file_edit(self):
        path = self.relation_file('mutable', ['SB(ab)=0'])
        public.configure_basis(front_rest=path)
        self.assertEqual(fb.NF_rels[2], 1)
        self.relation_file('mutable', ['SB(ab)=0', 'SB(ac)=0'], label='changedrels')
        self.assertEqual(fb.NF_rels[2], 1)
        self.assertEqual(public.relation_count(2), 1)
        public.clear_basis_cache()
        self.assertEqual(fb.NF_rels[2], 2)
        self.assertEqual(set(public.fr_rels(2)), {'ab', 'ac'})

    def test_weight_seven_uses_declared_block_before_legacy_fallback(self):
        path = self.relation_file('weight7', ['SB(aaaaaaa)=0'], weight=7)
        with path.open('a') as stream:
            stream.write('iheptrels_z := [SB(bbbbbbb)=0]:\n'
                         'iheptrels_c := [SB(ccccccc)=SB(ddddddd)]:\n')
        public.configure_basis(front_rest=path)
        self.assertEqual(fb.NF_rels[7], 1)
        self.assertEqual(public.fr_rels(7), {'aaaaaaa': {None: 0}})

    def test_weight_seven_still_supports_missing_combined_legacy_block(self):
        path = self.relation_file('legacy7', [], weight=7)
        path.write_text(path.read_text().replace('customrels := []:\n', '') +
                        'iheptrels_z := [SB(bbbbbbb)=0]:\n'
                        'iheptrels_c := [SB(ccccccc)=SB(ddddddd)]:\n')
        public.configure_basis(front_rest=path)
        self.assertEqual(fb.NF_rels[7], 2)
        self.assertEqual(set(public.fr_rels(7)), {'bbbbbbb', 'ccccccc'})

    def test_weight_zero_and_weights_beyond_old_tables(self):
        self.assertEqual(fb.NF_rels[0], 0)
        self.assertEqual(fb.NB_rels[0], 0)
        self.assertEqual(public.fr_rels(0), {})
        path = self.relation_file('weight12', ['SB(aaaaaaaaaaaa)=0'], weight=12)
        public.configure_basis(front_rest=path)
        self.assertEqual(fb.NF_rels[12], 1)

    def test_failed_load_is_not_cached_as_zero(self):
        path = self.root / 'missing'
        with self.assertRaises(FileNotFoundError):
            public.relation_count(2, prefix=path)
        self.assertFalse(fb.relation_counts)
        self.relation_file('missing', ['SB(ab)=0'])
        with self.assertRaises(KeyError):
            public.relation_count(3, prefix=path)
        self.assertFalse(fb.relation_counts)
        self.assertEqual(public.relation_count(2, prefix=path), 1)

    def test_invalid_queries_are_rejected(self):
        for weight in (-1, 1.5, True):
            with self.subTest(weight=weight):
                with self.assertRaises(ValueError):
                    public.relation_count(weight)
                with self.assertRaises(ValueError):
                    fb.NF_rels[weight]
        with self.assertRaises(ValueError):
            public.relation_count(2, 'unknown')
        with self.assertRaises(ValueError):
            public.relation_count(2, 'back', optrace='unknown')


if __name__ == '__main__':
    unittest.main()
