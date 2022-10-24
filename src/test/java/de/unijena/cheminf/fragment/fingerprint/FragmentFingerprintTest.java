/*
 * MIT License
 *
 * Copyright (c) 2022
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

package de.unijena.cheminf.fragment.fingerprint;

import org.junit.Assert;
import org.junit.Test;
import java.util.ArrayList;
import java.util.HashMap;

/**
 *  Class to test the correct working of BitFragmentFingerprint and CountFragmentFingerprint
 */
public class FragmentFingerprintTest extends BitFragmentFingerprint {

    //<editor-fold desc="Constructor">
    /**
     * Empty Constructor
     */
    public FragmentFingerprintTest() {
    }
    //</editor-fold>

    /**
     * Tests the instantiation and the call of the method "generateFingerprint"
     *
     * @throws Exception  if anything goes wrong
     */
    @Test
    public void basicTest() throws Exception {
        ArrayList<String> tmpTestList = new ArrayList<>();
        tmpTestList.add("*OS(=O)(=O)O[H]");
        tmpTestList.add("CC1(C)CCCC2(C)C1CCC3(C)C2CCC4C5CCCC5(C)CCC43C");
        tmpTestList.add("*N(*)S(*)(=O)=O");
        tmpTestList.add("*C(=O)C=C");
        tmpTestList.add("*OCO[H]");
        tmpTestList.add("CC(C)C1CCC2(C)C1CCC3(C)C2CCC4C5(C)CCCC(C)(C)C5CCC43C");
        tmpTestList.add("ccc(cCC(CC)C1CCC(CC)CC1)CC(C)CC2CC3(CCCC3C)CCC(C)(CC4CCCC(C4)C(C)C5CCCCC5)C2C(C)C6CCCC(Cc7ccccc7)C6");
        tmpTestList.add("c1cc(c(cc1CCC(C)CCCC(C)CCCC(C)CCCC(C)C)C)C");
        tmpTestList.add("*C(=O)N(*)C=CC(=O)N(*)*");
        tmpTestList.add("c1ccc(cc1)CCc2ccccc2CC");
        tmpTestList.add("O=C");
        tmpTestList.add("[H]OC1OC(=O)C=C1");
        tmpTestList.add("c1cc(cc(c1)CC(C)C(C)C(C)C)CCCCC2CC(CC)C2c3cccc(c3)CCC");
        tmpTestList.add("*OC=C(O[H])C(*)=O");
        tmpTestList.add("c1cc(cc(c1)CC)C");
        tmpTestList.add("CCCCCCCCCCCCCCCCCCCCCCC");
        tmpTestList.add("*C(=O)OC=CO[H]");
        tmpTestList.add("[H]OC=NC=C");
        tmpTestList.add("*OC=CC(=O)O*");
        tmpTestList.add("*N(*)*");
        tmpTestList.add("*OC(=CBr)C(Br)=C");
        tmpTestList.add("cc(c)C");
        tmpTestList.add("*s*");
        tmpTestList.add("cc");
        tmpTestList.add("CCCC1(C)C(C)CC(C)C(C)C1CC(CC)C(C)CC");
        tmpTestList.add("*OC=C(O*)C(*)=O");
        tmpTestList.add("c1ccc(cc1)CCC");
        tmpTestList.add("*OCN(*)C(=O)N=CN(*)*");
        tmpTestList.add("cC");
        tmpTestList.add("BrC=C");
        tmpTestList.add("*OP(=O)(O*)O[H]");


        BitFragmentFingerprint tmpBitFingerprintTest = new BitFragmentFingerprint();
        CountFragmentFingerprint tmpCountFingerprintTest = new CountFragmentFingerprint();
        HashMap<String,String> tmpMoleculeFragments = new HashMap<>();

        tmpMoleculeFragments.put("BrC=C", "25");
        tmpMoleculeFragments.put("ccc(cCC(CC)C1CCC(CC)CC1)CC(C)CC2CC3(CCCC3C)CCC(C)(CC4CCCC(C4)C(C)C5CCCCC5)C2C(C)C6CCCC(Cc7ccccc7)C6", "25");
        tmpMoleculeFragments.put("*s*", "10");
        tmpMoleculeFragments.put("*OC(C(*)=O)N(*)C(*)=O", "5");
        tmpMoleculeFragments.put("c1cc(cc(c1)C)c2cccc(c2)C3C4C(C)CCC4C3(C)C", "7");
        tmpMoleculeFragments.put("c1cc(ccc1C)C", "3");
        tmpMoleculeFragments.put("*OC(=O)C=CN(*)C(*)=O", "66");
        tmpMoleculeFragments.put("CCCC(C)C", "74");


        ArrayList tmpBitTest1 =  tmpBitFingerprintTest.generateFragmentFingerprint(tmpTestList, tmpMoleculeFragments);
        int[] tmpBiTVector1 = (int[]) tmpBitTest1.get(0);
        int[] tmpExceptedBitVector = {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0};
        Assert.assertArrayEquals(tmpExceptedBitVector,tmpBiTVector1);

        ArrayList tmpCountTest1 = tmpCountFingerprintTest.generateFragmentFingerprint(tmpTestList, tmpMoleculeFragments);
        int[] tmpCountVector1 = (int[]) tmpCountTest1.get(0);
        int[] tmpExceptedCountVector = {0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,0,25,0,0,0,0,0,0,0,0,0,0,25,0,0,0};
        Assert.assertArrayEquals(tmpExceptedCountVector,tmpCountVector1);
    }

}