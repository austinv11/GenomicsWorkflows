#!/usr/bin/env sh

# List /ref directory tree if no arguments are provided, otherwise run the given command
if [ $# -eq 0 ]; then
    tree /ref
else
    bash -c ". /bin/reference_functions.sh && export PATH=\"$PATH\"; export TENX_DISABLE_TELEMETRY=\"$TENX_DISABLE_TELEMETRY\"; $@ "
fi
