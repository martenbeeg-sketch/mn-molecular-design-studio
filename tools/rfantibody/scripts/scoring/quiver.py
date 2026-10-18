#!/usr/bin/env python3
'''
    This file defines the Quiver file class which is used to store PDB files and their associated scores.
    This class is going to be the simplest and quickest implementation of a database for the PDB files and their scores.
    This Quiver implementation will be just a list of PDB lines in a single file, with a tag for each PDB file.

    Later this can be made more sophisticated by using a proper database, but for now this will be the simplest implementation.
'''

import sys
import os

class Quiver():
    def __init__(self, filename, mode, backend='txt'):
        '''
            filename:   the name of the Quiver file to operate on
            mode:       the mode to open the file in, either 'r' for read-only, or 'w' for write-only
        '''

        self.mode = mode
        self.fn = filename

        self.backend = backend

        self.buffer = {}

        # Perform mode-specific operations
        if self.mode == 'w':
            self.tags = self._read_tags()
        elif self.mode == 'r':
            self.tags = self._read_tags()
        else:
            sys.exit(f'Quiver file must be opened in either read or write mode, not {self.mode}')

    def _read_tags(self) -> list:
        '''
            Read the tags from the Quiver file
        '''
        tags = []

        # Check if the file exists
        if not os.path.exists(self.fn):
            return tags

        with open(self.fn, 'r') as f:
            for line in f:

                if line.startswith('QV_TAG'):
                    tags.append(''.join(line.split()[1:]))

        return tags

    def get_tags(self):
        return self.tags

    def size(self):
        return len(self.tags)

    def add_pdb(self, pdb_lines: list, tag: str, score_str=None) -> None:
        '''
            Add a PDB file to the Quiver file
            The Quiver file must be opened in write mode to allow for writing.

            Inputs:
                pdb_lines:  a list of strings, each string is a line of the PDB file
                tag:        a string, the tag to associate with this PDB file
                score_str:  a string, the score to associate with this PDB file
        '''

        if self.mode == 'r':
            # We could in the future have this fail and return False
            sys.exit(f'Quiver file must be opened in write mode to allow for writing.')

        if tag in self.tags:
            # This can be made more sophisticated later
            sys.exit(f'Tag {tag} already exists in this file.')

        # We could eventually have a buffering system here to avoid writing to disk every time
        with open(self.fn, 'a') as f:
            f.write(f'QV_TAG {tag}\n')
            if score_str is not None:
                f.write(f'QV_SCORE {tag} {score_str}\n')
            f.write(''.join(pdb_lines))
            f.write('\n')
        
        self.tags.append(tag)

    def get_pdblines(self, tag: str) -> list:
        '''
            Get the PDB lines associated with the given tag
            The Quiver file must be opened in read mode to allow for reading.
            This function will iterate through the file until it finds the tag, then it will return the PDB lines associated with that tag.

            Inputs:
                tag:    a string, the tag to get the PDB lines for

            Outputs:
                pdblines:  a list of strings, each string is a line of the PDB file
        '''

        if self.mode == 'w':
            # We could in the future have this fail and return False
            sys.exit(f'Quiver file must be opened in read mode to allow for reading.')
        
        with open(self.fn, 'r') as f:
            for line in f:
                if line.startswith('QV_TAG'):
                    if tag == line.split()[1]:
                        pdb_lines = []
                        for line in f:
                            if line.startswith('QV_SCORE'):
                                continue
                            if line.startswith('QV_TAG'):
                                break
                            pdb_lines.append(line)

                        return pdb_lines

        # If we get here, we didn't find the tag
        sys.exit(f'Requested tag: {tag} which does not exist')

