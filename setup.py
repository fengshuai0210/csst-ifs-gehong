import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="csst-ifs-gehong",                                     
    version="0.0.1",                                        
    author="Shuai Feng",                                    
    author_email="sfeng@hebtu.edu.cn",                      
    description="A Python package for IFS data modelling",                            
    long_description=long_description,                      
    long_description_content_type="text/markdown",          
    url="https://github.com/",                              
    packages=setuptools.find_packages(),                    
    classifiers=[                                           
        "Programming Language :: Python :: 3",              
        "License :: OSI Approved :: MIT License",           
        "Operating System :: OS Independent",               
    ],
)