Installation
============

.. _installation:

Installation
------------

To use ``csst-ifs-gehong``, first install it using ``pip``:

.. code-block:: console

    $ pip install csst-ifs-gehong

Environment variable setting
----------------------------

The modelling of spectrum need some dataset such as spectral templates, 
which are not included in the `pip` installed package. You need to download 
those data separately from https://pan.baidu.com/s/1zg9IbxZU9sc63KvJQ1nWtg with
`wx79`. Unzip the `gehong_data.zip` file after downloading. Then, add the path where 
the `data` folder resides in the `gehong_data` to the environment variable `GEHONG_DATA_PATH`. 

For example, you downloaded the ``gehong_data.zip`` to ``~/GIT/``, and ``unzip`` it. 

.. code-block:: console

    sfeng_MacBook:GIT sfeng$ ls
    gehong_data.zip
    sfeng_MacBook:GIT sfeng$ unzip gehong_data.zip 
    ...
    sfeng_MacBook:GIT sfeng$ ls
    __MACOSX        gehong_data     gehong_data.zip

