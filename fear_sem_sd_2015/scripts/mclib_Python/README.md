# McIntyre Library
This is a set of simple classes to handle various types of 
genomic/transcriptomic data. Some of these classes are wrappers for 
other python modules such as `gffutils`. In order to use these classes various 
python modules will need to be installed. Including:

* cPickle
* numpy
* pandas
* matplotlib
* pysam
* gffutils
* pyvcf

I suggest using this repository as a subtree in your project. This will allow
you to pull changes to mclib_Python while the project is active. But then have
a code freeze once the project is published, so that you can easily tell what
code was run. First setup you project's scripts folder as a git repository.
Then to subtree this library do the following inside of your scripts folder.

**Add the mclib_Python remote**
```bash
$ git remote add mclib_Python git@github.com:McIntyre-Lab/mclib_Python.git
```

**Create the Subtree**
```bash
$ git subtree add --prefix=mclib_Python --squash mclib_Python master
```

**To pull updates from upstream**
```bash
$ git subtree pull --prefix=mclib_Python --squash mclib_Python master
```

**To push updates to upstream**
```bash
$ git subtree push --prefix=mclib_Python --squash mclib_Python master
```

For more detailed information about this library please see [McIntyre Library
Documentation](http://bio.rc.ufl.edu/pub/mcintyre/mcpython/mclib/index.html)
