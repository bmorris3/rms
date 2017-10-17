from setuptools import setup

setup(name='rms',
      version='0.1',
      description='Rotational Modulation Saturation calculator',
      install_requires=['numpy', 'astropy', 'matplotlib', 'scipy',
                        'sphinx-automodapi'],
      author='Brett Morris',
      author_email='bmmorris@uw.edu',
      license='MIT',
      url='https://github.com/bmorris3/mrspoc',
      zip_safe=False,
      use_2to3=False,
)
