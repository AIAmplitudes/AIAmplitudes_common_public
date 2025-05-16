

To package the repo:

git tag -a {version} -m "tag message"
hatch build
hatch publish

From AIAA_common_dev, do:

docker build
