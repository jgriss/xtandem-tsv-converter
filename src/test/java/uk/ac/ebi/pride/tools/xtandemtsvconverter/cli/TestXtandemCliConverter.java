package uk.ac.ebi.pride.tools.xtandemtsvconverter.cli;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.tools.xtandemtsvconverter.util.PSM;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created by jg on 16.07.15.
 */
public class TestXtandemCliConverter {
    File testFile;

    @Before
    public void setUp() throws Exception {
        testFile = new File(TestXtandemCliConverter.class.getClassLoader().getResource("PRD000001.mgf.xtandem.2015_07_15_21_15_03.t.xml").getFile());

        //testFile = new File("/tmp/classical_engine.xml");
    }

    @Test
    public void testGetPsms() throws Exception {
        List<PSM> psms = XtandemCliConverter.extractPsms(testFile);
        psms = new ArrayList<PSM>(psms);
        Collections.sort(psms);

        final double fdr = 0.01;

        int target = 0, decoy = 0;
        int aboveFdr = 0;

        for (PSM psm : psms) {
            if (psm.isDecoy())
                decoy++;
            else
                target++;

            double currentFdr = (double) decoy / target;

            if (currentFdr < fdr)
                aboveFdr++;
        }

        Assert.assertEquals(185, psms.size());

        PSM psm10 = psms.get(10);
        Assert.assertTrue(psm10.getSpectrumIndex() == 5746 || psm10.getSpectrumIndex() == 2775); // both identified as same with same expect
        Assert.assertEquals(0.013, psm10.getExpect(), 0.00000000001);
        Assert.assertEquals("CISHEYR", psm10.getSequence());
        Assert.assertFalse(psm10.isDecoy());
    }
}
