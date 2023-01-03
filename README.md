### Notices

#### Mirrors

Repository:
- [Codeberg](https://codeberg.org/paveloom-university/Resonances-in-Celestial-Mechanics-S11-2022)
- [GitHub](https://github.com/paveloom-university/Resonances-in-Celestial-Mechanics-S11-2022)
- [GitLab](https://gitlab.com/paveloom-g/university/s11-2022/resonances-in-celestial-mechanics)

#### Prerequisites

Make sure you have installed

- [Julia](https://julialang.org)
- TeX Live:
    - Packages (Fedora Linux):
        - `poppler-utils`
        - `texlive-pgfplots`
        - `texlive-standalone`
        - `texlive-xetex`

#### Run

First, instantiate the project with

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

Get the data by running the `data.bash` script.

Then, run one of the following

```bash
# Run without a daemon
julia --project=. scripts/script.jl
julia --project=. scripts/script.jl --postfix "'Custom postfix'"

# Run with a daemon
./julia.bash scripts/script.jl
./julia.bash scripts/script.jl --postfix "'Custom postfix'"
```

Use the `--help` flag to get help information.

To kill the daemon, run `./julia.bash kill`.
