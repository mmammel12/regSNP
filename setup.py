import setuptools

setuptools.setup(
    name="regsnp",
    version="0.2.8",
    packages=setuptools.find_packages(),
    scripts=["regsnp_intron/regsnp_intron.py"],
    install_requires=["pandas==0.17.1", "pysam==0.15.2", "pymongo"],
    include_package_data=True,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    url="https://github.com/mmammel12/regSNP",
    license="MIT",
    author="linhai, mamammel",
    author_email="linhai@iupui.edu, mamammel@iu.edu",
    description="Predict disease-causing probability of human intronic SNVs.",
)
