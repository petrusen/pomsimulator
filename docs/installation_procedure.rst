Installation procedure
======================

The system administrator must perform following actions prior deploying the ioChem-BD platform.

Create iochembd system user
---------------------------

First we will create a user account called *iochembd*, it will hold ioChem-BD files and will be responsible for running the web service that provides access to the software.

.. code:: bash

     root# useradd -m iochembd

We can use any other non privileged system user or username to install or run this software, we called it *iochembd* as an example.

.. attention:: The following commands will be run as the *iochembd* user, unless otherwise mentioned.