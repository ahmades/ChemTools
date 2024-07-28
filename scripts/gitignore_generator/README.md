# gitignore generator

## Purpose

Concatenates a comma-separated list of .gitignore configurations from https://github.com/github/gitignore and moves it to a specific destination, usually the root directory of a repository.

## Usage

If ran from the root directory of the repository, to include the .gitignore of `C++`, `Python`, `Visual Studio`, and `Visual Studio Code`, use:

```bash
./scripts/gitignore_generator/generate_gitignore.sh \
  --files C++.gitignore,Python.gitignore,VisualStudio.gitignore,Global/VisualStudioCode.gitignore \
  --destination . \
  --work_dir ~/work/tmp/gitignore
```

All arguments are mandatory.
