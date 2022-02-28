import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="tagee",
    version="0.1.0",
    description="A suite of terrain analysis tools for Google Earth Engine.",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/zecojls/tagee",
    author="Jose Lucas",
    author_email="zecojls@gmail.com",
    license_files=["LICENSE"],
    packages=["tagee"],
    include_package_data=True,
    install_requires=["earthengine-api"],
)