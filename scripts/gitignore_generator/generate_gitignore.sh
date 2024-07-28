#!/bin/bash

# Parse the arguments
parse_args() {
    while [[ $# -gt 0 ]]; do
        key="$1"
        case ${key} in
        --files)
            IFS=',' read -r -a FILES <<<"$2"
            shift 2
            ;;
        --destination)
            DESTINATION=$(realpath "$2")
            shift 2
            ;;
        --work_dir)
            WORK_DIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: ${key}"
            exit 1
            ;;
        esac
    done

    # Ensure all required arguments are provided
    if [ -z "$WORK_DIR" ] || [ -z "$DESTINATION" ] || [ ${#FILES[@]} -eq 0 ]; then
        echo "Usage: $0 --work_dir <work_dir> --destination <destination_path> --files <file1,file2,...>"
        exit 1
    fi
}

# Creates work directory
function create_work_dir() {
    local work_dir=$1
    mkdir -p "{$work_dir}"
}

# Removes work directory
remove_work_dir() {
    local work_dir=$1
    rm -rf "$work_dir"
}

# Clones the gitignore repository to the work directory
clone_repo() {
    local work_dir=$1
    git clone https://github.com/github/gitignore "${work_dir}/gitignore"
    if [ $? -ne 0 ]; then
        echo "Failed to clone the repository."
        exit 1
    fi
}

# Concatenate any number of .gitignore templates into a single .gitignore file
concatenate_gitignore() {
    local work_dir=$1
    shift
    local gitignore_files=("$@")
    local output_file="$work_dir/.gitignore"
    >"$output_file"
    for file in "${gitignore_files[@]}"; do
        local file_path="$work_dir/gitignore/$file"
        if [ ! -f "$file_path" ]; then
            echo "Required .gitignore file $file_path not found."
            exit 1
        fi
        echo "Adding $file"
        echo "###############################################" >>"$output_file"
        echo "# $file" >>"$output_file"
        echo "# source: https://github.com/github/gitignore" >>"$output_file"
        echo "###############################################" >>"$output_file"
        cat "$file_path" >>"$output_file"
        echo -e "\n" >>"$output_file"
    done
}

# Moves the generated .gitignore file to the user-defined destination path
move_gitignore() {
    local work_dir=$1
    local destination=$2
    mv "${work_dir}/.gitignore" "${destination}"
    if [ $? -ne 0 ]; then
        echo "Failed to move the .gitignore file to the destination."
        exit 1
    fi
}

main() {
    parse_args "$@"
    create_work_dir "${WORK_DIR}"
    clone_repo "${WORK_DIR}"
    concatenate_gitignore "${WORK_DIR}" "${FILES[@]}"
    move_gitignore "${WORK_DIR}" "${DESTINATION}"
    remove_work_dir "${WORK_DIR}"
    echo "Successfully created and moved .gitignore to ${DESTINATION}"
}

main "$@"
