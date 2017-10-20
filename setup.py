from setuptools import setup
import os
import sys
from setuptools.command.install import install
import subprocess

packagename = 'rms'


def get_virtualenv_path():
    """Used to work out path to install compiled binaries to."""
    if hasattr(sys, 'real_prefix'):
        return sys.prefix

    if hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix:
        return sys.prefix

    if 'conda' in sys.prefix:
        return sys.prefix

    return None


def compile_and_install_software():
    """
    Used the subprocess module to compile/install the C software.

    Return the path of the new STSP executable, and the path we want to
    move it to when the installation is complete.
    """
    src_path = './STSP/'

    # compile the software
    cmd = "gcc -lm stsp.c -o stsp_rms"
    venv = get_virtualenv_path()
    if venv:
        cmd += ' --prefix=' + os.path.abspath(venv)
    subprocess.check_call(cmd, cwd=src_path, shell=True)
    pkg_dir = os.path.join(sys.exec_prefix, 'lib',
                           "python{0}.{1}".format(sys.version_info.major,
                                                  sys.version_info.minor),
                           'site-packages', packagename)
    executable_build_path = os.path.join(src_path, 'stsp_rms')
    executable_new_path = os.path.join(pkg_dir, 'stsp_rms')
    return executable_build_path, executable_new_path


class CustomInstall(install):
    """Custom handler for the 'install' command."""
    def run(self):
        #paths = compile_and_install_software()
        super().run()
        #os.rename(*paths)


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
      cmdclass={'install': CustomInstall}
)
