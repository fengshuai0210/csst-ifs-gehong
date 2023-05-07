from setuptools import setup, find_packages


setup(
    name='csst-ifs-gehong',
    version='1.0.0',
    license='MIT',
    author="Shuai Feng",
    author_email='sfeng@hebtu.edu.cn',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    url='https://csst-ifs-gehong.readthedocs.io/en/latest/',
    keywords='CSST-IFS',
    install_requires=[
          'astropy',
      ],

)