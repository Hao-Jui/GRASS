# CI parameter files

CTest and the Make test targets compile the solver against
`tests/CI/para_panel.f90`. Production executables continue to compile against
`src/para_panel.f90`, so changing run parameters cannot silently change a
regression fixture.

Keep the public `para_mod` schema and allocation routines identical in both
files. Change values in the CI copy only when intentionally updating a test
case and its reference output.
