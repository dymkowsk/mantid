from __future__ import (absolute_import, division, print_function)


import unittest
from mantid.simpleapi import *
from mantid.api import *


class IndirectQuickDiffractionTest(unittest.TestCase):

    def test_basic(self):
        """
        Sanity test to ensure the most basic reduction actually completes.
        """

        IndirectQuickReduction(RunNumbers='26176', Instrument='IRIS', Analyser='graphite', Reflection='002', SpectraRange=[3, 53],
                               ElasticRange=[-0.5, 0], InelasticRange=[0, 0.5], GroupingMethod='All', DiffractionSpectra=[105, 112])

        wks = mtd['iris26716_diffspec']


        self.assertTrue(isinstance(wks, WorkspaceGroup), 'Result workspace should be a workspace group.')
        self.assertEqual(len(wks), 1)
        self.assertEqual(wks.getNames()[0], 'iris26176_diffspec_red')

        red_ws = wks[0]
        self.assertEqual(red_ws.getAxis(0).getUnit().unitID(), 'dSpacing')
        self.assertEqual(red_ws.getNumberHistograms(), 1)
