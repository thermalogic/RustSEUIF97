
from setuptools import find_packages, setup
from setuptools_rust import RustExtension

setup(
    packages=find_packages(where="python"),
    package_dir={"": "python"},
    rust_extensions=[RustExtension("seuif97.seuif97",
                                   features=["python"],)],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
    ]
)
