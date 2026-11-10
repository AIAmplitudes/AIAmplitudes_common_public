"""Data-driven basis dimensions, source isolation, and cache lifecycle."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import aiamplitudes_common_public as public
from aiamplitudes_common_public.fbspaces import (
    F_number, B_number, F_perm_number, B_perm_number,
    basis_sizes, clear_basis_cache,
)


class BasisSizeTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.sources = (public.frontspace_file, public.frontspace_perm_file,
                        public.backspace_perm_file, dict(public._backspace_files))
        clear_basis_cache()

    def tearDown(self):
        (public.frontspace_file, public.frontspace_perm_file,
         public.backspace_perm_file, backfiles) = self.sources
        public._backspace_files.update(backfiles)
        clear_basis_cache()

    def restrictive_file(self, filename, words, weight=2, space='front'):
        path = self.root / filename
        header, wrapper = ('init', 'SB') if space == 'front' else ('fin', 'E')
        entries = ','.join(f'{wrapper}({word})' for word in words)
        path.write_text(f'{header}_list[{weight}] := testindep :\n'
                        f'testindep := [{entries}]:\n')
        return path

    def permissive_file(self, filename, words, space='front'):
        path = self.root / filename
        entries = ','.join(f'SB({word})' for word in words)
        path.write_text(f'{space}space[2] := [{entries}]:\n')
        return path

    def test_load_and_size_lookup_share_one_load_and_one_count(self):
        class CountedBasis(dict):
            counts = 0

            def __len__(self):
                self.counts += 1
                return super().__len__()

        basis = CountedBasis({'FR_2_1': 'aa', 'FR_2_2': 'bb'})
        with patch.object(public, 'get_rest_fspace', return_value=(basis, {})) as loader:
            self.assertIs(public.fspace(2, 'R'), basis)
            for _ in range(5):
                self.assertEqual(public.basis_size(2), 2)
                self.assertEqual(F_number[2], 2)
                self.assertIs(public.fspace(2, 'R'), basis)
                public.fspace_flip(2, 'R')
            self.assertEqual(loader.call_count, 1)
            self.assertEqual(basis.counts, 1)

    def test_files_at_the_same_weight_keep_separate_sizes(self):
        first = self.restrictive_file('first', ['aa', 'bb'])
        second = self.restrictive_file('second', ['aa', 'bb', 'cc', 'ab'])
        for path, size in ((first, 2), (second, 4), (first, 2)):
            self.assertEqual(len(public.fspace(2, 'R', prefix=path)), size)
            self.assertEqual(public.basis_size(2, prefix=path), size)
        self.assertEqual(len(basis_sizes), 2)

    def test_selected_sources_drive_legacy_views_and_keep_forms_separate(self):
        public.configure_basis(
            front_rest=self.restrictive_file('fr', ['aa', 'bb']),
            back_rest=self.restrictive_file('br', ['aa'], space='back'),
            front_perm=self.permissive_file('fp', ['aa', 'bb', 'cc']),
            back_perm=self.permissive_file('bp', ['aa', 'bb'], space='back'),
        )
        self.assertEqual((F_number[2], B_number[2], F_perm_number[2], B_perm_number[2]),
                         (2, 1, 3, 2))
        public.configure_basis(front_rest=self.restrictive_file('new', ['cc']))
        self.assertEqual(F_number[2], 1)
        self.assertEqual(public.fspace(2, 'R'), {'FR_2_1': 'cc'})
        self.assertEqual(F_perm_number[2], 3)

    def test_operator_traces_do_not_share_dimensions(self):
        public.configure_basis(back_rest=self.restrictive_file('phi2', ['aa'], space='back'))
        public.configure_basis(back_rest=self.restrictive_file('phi3', ['aa', 'bb'], space='back'),
                               optrace='phi3')
        self.assertEqual(public.basis_size(2, 'back', optrace='phi2'), 1)
        self.assertEqual(public.basis_size(2, 'back', optrace='phi3'), 2)
        self.assertEqual(B_number[2], 1)

    def test_weights_beyond_old_tables_and_weight_zero(self):
        path = self.restrictive_file('weight20', ['a' * 20, 'b' * 20], weight=20)
        public.configure_basis(front_rest=path)
        self.assertEqual(F_number[20], 2)
        self.assertEqual(F_number[0], 1)
        self.assertEqual(public.fspace(0, 'R'), {'FR_0_1': ''})
        self.assertEqual(public.basis_size(0, 'back', 'P'), 1)

    def test_explicit_clear_reloads_file_contents_and_header(self):
        path = self.restrictive_file('mutable', ['aa'])
        public.configure_basis(front_rest=path)
        self.assertEqual(F_number[2], 1)
        self.restrictive_file('mutable', ['aa', 'bb', 'cc'])
        self.assertEqual(F_number[2], 1)
        clear_basis_cache()
        self.assertEqual(F_number[2], 3)

    def test_weight_six_special_loader_and_flip_share_cached_basis(self):
        basis = {'FP_6_1': {'aaaaaa': 1}, 'FP_6_2': {'bbbbbb': 2}}
        with patch.object(public, 'get_perm_fspace_wt6', return_value=basis) as loader:
            self.assertIs(public.fspace(6, 'P'), basis)
            self.assertEqual(F_perm_number[6], 2)
            self.assertEqual(public.fspace_flip(6, 'P'), {
                'aaaaaa': {'FP_6_1': 1}, 'bbbbbb': {'FP_6_2': 2},
            })
            self.assertEqual(loader.call_count, 1)

    def test_failed_load_does_not_cache_a_size(self):
        with self.assertRaises(AssertionError):
            public.basis_size(2, prefix=self.root / 'missing')
        self.assertFalse(basis_sizes)


if __name__ == '__main__':
    unittest.main()
