# Add current directory to PATH and LD_LIBRARY_PATH in the user's shell rc file
# English comments only, code and comments must be in English

# Detect the current shell type
current_shell=$(basename "$SHELL")

# Set the shell rc file based on the shell type
if [ "$current_shell" = "bash" ]; then
    shell_rc="$HOME/.bashrc"
elif [ "$current_shell" = "zsh" ]; then
    shell_rc="$HOME/.zshrc"
elif [ "$current_shell" = "ksh" ]; then
    shell_rc="$HOME/.kshrc"
elif [ "$current_shell" = "csh" ]; then
    shell_rc="$HOME/.cshrc"
elif [ "$current_shell" = "tcsh" ]; then
    shell_rc="$HOME/.tcshrc"
else
    # Default to .profile if shell is unknown
    shell_rc="$HOME/.profile"
fi

# Get the absolute path of the current directory
current_dir="$(pwd)"

# Prepare export lines
export_path_line="export PATH=\"\$PATH:$current_dir\""
export_ld_line="export LD_LIBRARY_PATH=\"\$LD_LIBRARY_PATH:$current_dir/lib\""

# Append to shell rc file if not already present
if ! grep -Fxq "$export_path_line" "$shell_rc"; then
    echo "$export_path_line" >> "$shell_rc"
fi

if ! grep -Fxq "$export_ld_line" "$shell_rc"; then
    echo "$export_ld_line" >> "$shell_rc"
fi

# Source the shell rc file to apply changes immediately
# Use the correct command for the shell
if [ "$current_shell" = "csh" ] || [ "$current_shell" = "tcsh" ]; then
    source_cmd="source $shell_rc"
else
    source_cmd=". $shell_rc"
fi

eval "$source_cmd"
