from setuptools import setup,find_packages

with open("README.rst", "r", encoding="utf-8") as fh:
        long_description = fh.read()

setup(name='pomsimulator',
      version='1.0.1',
      author="Enric Petrus, Jordi Buils. Diego Garay-Ruiz",
      author_email="enricpz@icloud.com, jbuils@iciq.es, dgaray@iciq.es",
      description="Simulate the aqueous speciation of polyoxometalates (POMs) from quantum mechanical results",
      packages=find_packages(),
      install_requires=['numpy>=1.17.3','matplotlib>=3.1.2','networkx>=2.4','scipy>=1.6.1','pandas>=1.5', 'scikit-learn>=1.3.2','seaborn'],
      long_description=long_description,
      long_description_content_type="text/x-rst",  # Specify the content type as reStructuredText
      url="https://github.com/petrusen/pomsimulator",
      classifiers=[ "Programming Language :: Python :: 3",
                    "License :: OSI Approved :: MIT License",
                    "Operating System :: OS Independent"],
      python_requires='>=3.8',
      )
