#!/usr/bin/env bash
DIR="$(cd "$(dirname "$0")" && pwd)"
python -c "import sys; assert sys.version_info >= (3, 6), 'Sorry, you need python 3.6. If youre on midway, there is a conda env at $PRODUCTS/miniconda for you to use'" && python $DIR/submit_batch/submit_batch_jobs.py "$@"
