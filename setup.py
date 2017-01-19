from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(
    name = "nanoraw",
    version = "0.4.1",
    packages = ["nanoraw"],
    install_requires = ['h5py', 'numpy', 'scipy'],

    author = "Marcus Stoiber",
    author_email = "mhstoiber@lbl.gov",
    description='Analysis of nanopore sequencing data.',
    long_description = readme(),
    license = "BSD",
    keywords = "nanopore high-throughput sequencing correction genome",
    url = "https://github.com/marcus1487/nanoraw",

    entry_points={
        'console_scripts': [
            'nanoraw = nanoraw.__main__:main'
        ]
    },
    include_package_data=True,

    test_suite='nose2.collector.collector',
    tests_require=['nose2'],
)
