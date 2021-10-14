#!/usr/bin/env bash
DIR="$(cd "$(dirname "$0")" && pwd)"
python -c "import sys; assert sys.version_info >= (3, 6), 'Sorry, you need python 3.6. If youre on midway'" && python $DIR/makeDataFiles/makeDataFiles_main.py "$@"
