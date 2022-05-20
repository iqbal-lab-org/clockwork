# Clockwork multipass Virtual Machine

This is for developers to create a VM to develop the code and run the
tests. Assumes you have multipass installed.

To create a new machine called `VM_NAME` that has all the dependencies:

```
./launch_vm.sh VM_NAME
```

Log into the machine with:

```
multipass shell VM_NAME
```

The VM has your host $HOME directory mounted as $HOME/Home in the VM.

You can run the tests by navigating to the `python/` directory and running:

```
python3 setup.py test
```
