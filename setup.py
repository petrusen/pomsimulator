from setuptools import setup,find_packages
from os import path

here = path.abspath(path.dirname(__file__))

with open("README.rst", "r", encoding="utf-8") as fh:
        long_description = fh.read()

with open(path.join(here, "requirements.txt")) as requirements_file:
        # Parse requirements.txt, ignoring any commented-out lines.
        requirements = [line for line in requirements_file.read().splitlines() if not line.startswith("#")]


setup(name='pomsimulator',
      version='1.0.1a8',
      author="Enric Petrus, Jordi Buils. Diego Garay-Ruiz",
      author_email="enricpz@icloud.com, jbuils@iciq.es, dgaray@iciq.es",
      description="Simulate the aqueous speciation of polyoxometalates (POMs) from quantum mechanical results",
      packages=find_packages(),
      install_requires=requirements,
      long_description=long_description,
      long_description_content_type="text/x-rst",  # Specify the content type as reStructuredText
      url="https://github.com/petrusen/pomsimulator",
      classifiers=[ "Programming Language :: Python :: 3",
                    "License :: OSI Approved :: MIT License",
                    "Operating System :: OS Independent"],
      python_requires='>=3.8',
      include_package_data=True,
      package_data={ "inputs": [  # When adding files here, remember to update MANIFEST.in as well,
                                  # or else they will not be included in the distribution on PyPI!
                                  # 'path/to/data_file',
                                  "W_Set_PBE_test/*out", "W_Set_PBE_molfiles/*mol"],
                     "utilities": ["np_IM_test.csv"]}
      )
