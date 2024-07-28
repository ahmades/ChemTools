#!/bin/bash

function catch_error() {
    local exit_code=$1
    if [ "${exit_code}" != "0" ]; then
        echo "Error occurred"
        echo "  Line     : ${BASH_LINENO[1]}"
        echo "  Function : ${FUNCNAME[1]}"
        echo "  Command  : ${BASH_COMMAND}"
        echo "  Exit code: ${exit_code}"
    fi
}

function parse_arguments() {
    while [[ $# -gt 0 ]]; do
        key="$1"
        case ${key} in
        --clang_format_file)
            declare -g clang_format_file
            clang_format_file=$(realpath "$2")
            shift 2
            ;;
        --directories)
            shift
            # Read all arguments following --directories into the array
            #declare -A directories
            directories=()
            while [[ $# -gt 0 && ! "${key}" =~ ^-- ]]; do
                directories+=("$(realpath "${key}")")
                shift
            done
            ;;
        -* | *)
            echo "Unknown argument: ${key}"
            exit 1
            ;;
        esac
    done

    # required parameters
    local do_exit=false
    if [[ -z ${clang_format_file} ]]; then
        echo "Missing parameter --clang_format_file"
        do_exit=true
    fi

    if [[ -z "${directories[*]}" ]]; then
        echo "Missing parameter --directories"
        do_exit=true
    fi

    if ${do_exit}; then
        exit 1
    fi
}

function check_path_exists() {
    local path="$1"
    if ! [ -e "${path}" ]; then
        echo "${path} does not exists"
        exit 1
    fi
}

function validate_arguments() {
    check_path_exists "${clang_format_file}"

    for directory in "${directories[@]}"; do
        check_path_exists "${directory}"
    done
}

function apply_clang_format() {
    # Find all C++ source files with extensions .cpp and .hpp in the specified directories
    local cpp_files
    cpp_files=$(
        find "${directories[@]}" -type f \( \
            -name "*.C" -o \
            -name "*.cpp" -o \
            -name "*.cxx" -o \
            -name "*.c++" -o \
            -name "*.cc" -o \
            -name "*.h" -o \
            -name "*.hpp" -o \
            -name "*.hxx" -o \
            -name "*.hh" \
            \)
    )

    # Loop over each file and apply clang format
    for cpp_file in $cpp_files; do
        clang-format \
            --verbose \
            --Werror \
            --style=file:"${clang_format_file}" \
            -i \
            --assume-filename="${cpp_file}" \
            "${cpp_file}"
    done
}

function main() {
    trap 'catch_error $?' EXIT
    parse_arguments "$@"
    validate_arguments
    apply_clang_format
}

main "$@"
