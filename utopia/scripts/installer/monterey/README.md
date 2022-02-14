# Few notes on MacOS montery with M1 chips

## 10.02.2022

- spack fails to install any dependency which requires gfortran
- macports creates an environment that makes utopia_fe with the `stk backend` fail to link
- `brew` seems to be the only package manager that will work at this particular time

Issues with trilinos installation (to be done from source)
- use `git apply trilinos_monterey.patch` within trilinos source tree to fix linking issues