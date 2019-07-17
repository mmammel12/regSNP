try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name="regsnp_intron",
    version="0.2.0",
    packages=["regsnp_intron", "regsnp_intron.utils"],
    scripts=["bin/regsnp_intron"],
    install_requires=["pandas", "pysam", "pymongo"],
    include_package_data=True,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    url="https://github.com/linhai86/regsnp_intron",
    license="MIT",
    author="linhai",
    author_email="linhai@iupui.edu",
    description="Predict disease-causing probability of human intronic SNVs.",
)
