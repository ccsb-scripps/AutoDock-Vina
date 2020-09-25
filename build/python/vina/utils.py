#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Vina - utils
#

import os


def check_file_writable(fnm):
    """Source: https://www.novixys.com/blog/python-check-file-can-read-write/"""
    if os.path.exists(fnm):
        # path exists
        if os.path.isfile(fnm): # is it a file or a dir?
            # also works when file is a link and the target is writable
            return os.access(fnm, os.W_OK)
        else:
            return False # path is a dir, so cannot write as a file
    # target does not exist, check perms on parent dir
    pdir = os.path.dirname(fnm)
    if not pdir: pdir = '.'
    # target is creatable if parent dir is writable
    return os.access(pdir, os.W_OK)
