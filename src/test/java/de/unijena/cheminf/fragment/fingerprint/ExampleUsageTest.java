/*
 * MIT License
 *
 * Copyright (c) 2023 Betuel Sevindik, Felix Baensch, Jonas Schaub, Christoph Steinbeck, and Achim Zielesny
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

import org.junit.jupiter.api.Test;
import org.openscience.cdk.fragment.ExhaustiveFragmenter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;

import java.io.FileInputStream;
import java.util.HashMap;

public class ExampleUsageTest {
    /**
     *
     */
    @Test
    public void exampleUsageTest() throws Exception {
        //make it "path\\to\\coconut_napthalene_substructure_search_result.sdf" in the example code
        IteratingSDFReader tmpSDFReader = new IteratingSDFReader(new FileInputStream("D:\\Project_Fragment_Fingerprints\\FragmentFingerprints_git\\src\\test\\resources\\de\\unijena\\cheminf\\fragment\\fingerprint\\coconut_napthalene_substructure_search_result.sdf"), SilentChemObjectBuilder.getInstance());
        //This fragmentation scheme simply breaks single non-ring bonds.
        ExhaustiveFragmenter tmpFragmenter = new ExhaustiveFragmenter();
        tmpFragmenter.setMinimumFragmentSize(1); //minimum fragment size 1
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique);
        HashMap<String, Integer> tmpFrequenciesMap = new HashMap<>(50, 0.75f);
        while (tmpSDFReader.hasNext()) {
            IAtomContainer tmpMolecule = tmpSDFReader.next();
            tmpFragmenter.generateFragments(tmpMolecule);
            IAtomContainer[] tmpFragments = tmpFragmenter.getFragmentsAsContainers();
            for (IAtomContainer tmpFragment : tmpFragments) {
                String tmpSmilesCode = tmpSmiGen.create(tmpFragment);
                if (tmpFrequenciesMap.containsKey(tmpSmilesCode)) {
                    tmpFrequenciesMap.put(tmpSmilesCode, tmpFrequenciesMap.get(tmpSmilesCode) + 1);
                } else {
                    tmpFrequenciesMap.put(tmpSmilesCode, 1);
                }
            }
        }
        for (String tmpFragmentSmilesCode : tmpFrequenciesMap.keySet()) {
            System.out.println(tmpFragmentSmilesCode + ": " + tmpFrequenciesMap.get(tmpFragmentSmilesCode));
        }
        //TODO: set up fingerprinter
        String tmpCNP0170165SmilesString = "OC=1C(O)=C(O)C=2C=CC=CC2C1O"; //naphthalene-1,2,3,4-tetrol
        String tmpCNP0437667SmilesString = "O=C(O)C1=CC=2C=CC=CC2C=C1O"; //3-hydroxy-2-naphthoic acid
        //TODO: Generate and print fingerprints for the two molecules
    }
}
