try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='regsnp_intron',
    version='0.1.6',
    packages=['regsnp_intron', 'regsnp_intron.rbp', 'regsnp_intron.utils', 'regsnp_intron.predictor',
              'regsnp_intron.junc_score', 'regsnp_intron.conservation', 'regsnp_intron.protein_feature'],
    scripts=['bin/regsnp_intron'],
    install_requires=['numpy',
                      'scipy',
                      'pandas',
                      'sklearn',
                      'pysam',
                      'bx-python',
                      'pybedtools'],
    include_package_data=True,
    classifiers=['Development Status :: 3 - Alpha',
                 'License :: OSI Approved :: MIT License',
                 'Programming Language :: Python :: 2.7',
                 'Topic :: Scientific/Engineering :: Bio-Informatics'],
    url='https://github.com/linhai86/regsnp_intron',
    license='MIT',
    author='linhai',
    author_email='linhai@iupui.edu',
    description='Predict disease-causing probability of human intronic SNVs.'
)
