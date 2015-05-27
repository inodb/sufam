# So U Found A Mutation? (SUFAM)

Found a mutation in one or more samples? Now you want to check if they are in
another sample. Unfortunately mutect, varscan or whatever other variant caller
is not calling them. Use SUFAM. The super sensitive validation caller that calls
everything on a given position.

## Installation

```
git clone http://github.com/inodb/sufam
cd sufam
python setup.py install
```

## Run
```
sufam --help
```
