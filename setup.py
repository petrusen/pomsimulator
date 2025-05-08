from setuptools import setup,find_packages
setup(name='pomsimulator',
      version='2.0.1',
      author="Enric Petrus, Jordi Buils, Diego Garay-Ruiz",
      author_email="enricpz@icloud.com, jbuils@iciq.es, dgaray@iciq.es",
      description="Simulate the aqueous speciation of polyoxometalates (POMs) from quantum mechanical results",
      packages=find_packages(),
      install_requires=['numpy>=1.17.3','matplotlib>=3.1.2','networkx>=2.4','scipy>=1.6.1','pandas>=1.5',
                        'scikit-learn>=1.3.2','seaborn', 'setuptools']
      )
