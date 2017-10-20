from setuptools import setup

packagename = 'rms'

setup(name=packagename,
      version='0.1',
      description='Rotational Modulation Simulator',
      install_requires=['numpy', 'astropy', 'matplotlib', 'scipy'],
      author='Brett Morris',
      author_email='bmmorris@uw.edu',
      license='MIT',
      url='https://github.com/bmorris3/rms',
      zip_safe=False,
      use_2to3=False,
      packages=[packagename],
      include_package_data=True,
      package_data={"": ["LICENSE"]},
)
