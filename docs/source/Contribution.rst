Contribution
+++++++++++++

All kinds of contributions are welcome:

    * Improvement of code with new features, bug fixes, and  bug reports
    * Improvement of documentation
    * Additional tests

If you want to contribute to code:

    1. Fork the latest main branch
    2. Create a dev environment: install dependencies and install tesspy in `develop mode <https://python-packaging-tutorial.readthedocs.io/en/latest/setup_py.html#develop-mode>`_
    3. Write failing tests
    4. Write new code
    5. Run tests and make sure they pass
    6. Update documentation
    7. Format code
    8. Pull Request

If you have any ideas or questions, feel free to open an issue.



1. Fork tesspy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This consists of three main steps:

    1. Forking tesspy
    2. Cloning your forked repo
    3. Sync your fork with the original tesspy repo


A fork is simply a copy of a repository. It allows to freely
experiment with changes without having an effect on the original project.
The easiest way to fork a repository (first step) is by using the web interface
on GitHub.org:

    1. Navigate to `tesspy repository <https://github.com/siavash-saki/tesspy>`_.
    2. Click **Fork** in the top-right corner of the page.
    3. Follow the procedure (select owner, add a description, etc.)
    4. Click **Create fork**.

Then you need to clone your forked repo by using the command line::

    $ git clone https://github.com/YOUR-USERNAME/tesspy tesspy

Finally, you need to sync your repository to
the upstream (original project) tesspy repository::

    $ cd tesspy
    $ git remote add upstream https://github.com/siavash-saki/tesspy.git

Verify the new upstream repository by::

    $ git remote -v

For more information about fork, you can visit the
`GitHub docs <https://docs.github.com/en/get-started/quickstart/fork-a-repo>`_.

2. Create a dev env
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, you need to create a clean conda environment by::

    $ conda create -c conda-forge -n tesspy_dev

Then, the tesspy dependencies must be installed as follows::

    $ conda activate tesspy_dev
    $ conda install -c conda-forge geopandas scipy h3-py osmnx hdbscan mercantile scikit-learn

Next, navigate to the directory, where your clone of tesspy is located.
Install tesspy in develop mode using::

    $ python setup.py develop

This way, you have installed tesspy in develop mode. You can now track your changes
in real-time in your dev environment and do not need to install the package each time
you change something in the code.

3. Write failing tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is a good practice to start writing tests even before writing any
code. All the tests should initially fail. Think about your desired
feature and write corresponding test cases.
All tests are in the ``tests`` directory. New tests must also be saved here.
For more info on Test-driven development, take a look
`here <https://en.wikipedia.org/wiki/Test-driven_development>`_.

4. Write new code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Then it is time to start modifying the code, writing new functions, building
new features, etc. You should write to the point when your initially created
test cases pass.

5. Run tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You should then run the test suit inside your own clone of the repository
using ``pytest``.


6. Update documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

tesspy documentation is in the folder ``docs``. It is written using
reStructuredText, `which is explained here <http://www.sphinx-doc.org/en/stable/rest.html#rst-primer>`_.

After adding any code or feature, please also add to the documentation. Make sure to
check the documentation to build correctly by rendering it using ``sphinx``.

7. Format code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

tesspy uses the `PEP8 <http://www.python.org/dev/peps/pep-0008/>`_ standard and
``black`` to ensure a consistent code format throughout the project.

So, before committing your changes, format the code using ``black``.

it is a good idea to integrate ``black`` into your IDE. For example, it is explained
`here <https://black.readthedocs.io/en/stable/integrations/editors.html#pycharm-intellij-idea>`_
how to use black with ``pycharm``.


8. Pull Request
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
When you're finished making your changes and have made sure everything
is working properly, you can submit a Pull Request. You can find more information
on PRs `here <https://help.github.com/articles/using-pull-requests/>`_.



