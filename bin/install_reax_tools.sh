#!/bin/bash
# Install script for ReaxTools - adds bin directory to PATH and LD_LIBRARY_PATH
# This script ensures both reax_tools binary and reax_plot.py script are available in PATH

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

# Get the absolute path of the bin directory (where this script is located)
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Prepare export lines - add bin directory to PATH
export_path_line="export PATH=\"\$PATH:$script_dir\""
export_ld_line="export LD_LIBRARY_PATH=\"\$LD_LIBRARY_PATH:$script_dir/lib\""

# Remove old reax_tools PATH entries to ensure clean installation
# This removes any previous reax_tools PATH entries to avoid conflicts
temp_file=$(mktemp)
if [ "$current_shell" = "csh" ] || [ "$current_shell" = "tcsh" ]; then
    # For csh/tcsh, filter out lines containing reax_tools
    grep -v "reax_tools" "$shell_rc" > "$temp_file" 2>/dev/null || cat "$shell_rc" > "$temp_file" 2>/dev/null
else
    # For bash/zsh/ksh, filter out lines containing reax_tools
    grep -v "reax_tools" "$shell_rc" > "$temp_file" 2>/dev/null || cat "$shell_rc" > "$temp_file" 2>/dev/null
fi
mv "$temp_file" "$shell_rc"

# Append new export lines to shell rc file if not already present
if ! grep -Fxq "$export_path_line" "$shell_rc"; then
    echo "" >> "$shell_rc"
    echo "# ReaxTools installation - added $(date)" >> "$shell_rc"
    echo "$export_path_line" >> "$shell_rc"
fi

if ! grep -Fxq "$export_ld_line" "$shell_rc"; then
    echo "$export_ld_line" >> "$shell_rc"
    echo "# End of ReaxTools installation" >> "$shell_rc"
fi

# Make scripts executable
chmod +x "$script_dir/reax_tools" 2>/dev/null || true
chmod +x "$script_dir/reax_plot.py" 2>/dev/null || true
chmod +x "$script_dir/reax_network.py" 2>/dev/null || true
chmod +x "$script_dir/reax_network_solver.py" 2>/dev/null || true

echo "ReaxTools installation completed!"
echo ""
echo "Installed components:"
echo "  - reax_tools (main analysis binary)"
echo "  - reax_plot.py (plotting script)"
echo "  - reax_network.py (network analysis)"
echo "  - reax_network_solver.py (network solver)"
echo ""
echo "Installation directory: $script_dir"
echo "Shell configuration updated: $shell_rc"
echo ""

# Source the shell rc file to apply changes immediately
if [ "$current_shell" = "csh" ] || [ "$current_shell" = "tcsh" ]; then
    source_cmd="source $shell_rc"
else
    source_cmd=". $shell_rc"
fi

echo "Applying changes to current session..."
eval "$source_cmd"

echo ""
echo "You can now use 'reax_tools' and 'reax_plot.py' from any directory."
echo "If the commands are not available, please restart your terminal or run:"
echo "  source $shell_rc"
