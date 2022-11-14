/*
 * MIT License
 *
 * Copyright (c) 2022 Betuel Sevindik, Felix Baensch, Jonas Schaub, Christoph Steinbeck, and Achim Zielesny
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
 */

package de.unijena.cheminf.fragment.fingerprint;

import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.ICountFingerprint;
import org.openscience.cdk.fingerprint.IFingerprinter;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Interface for implementing fragment fingerprints
 *
 * @author Betuel Sevindik
 */
public interface IFragmentFingerprinter extends IFingerprinter {

    /**
     * Method to generate bit fingerprint all by itself by comparing the unique SMILES
     *
     * @param aListOfUniqueSmiles  a list containing all the fragments produced by the fragmentation of a given molecule
     * @return IBitFingerprint
     */
   IBitFingerprint getBitFingerprint(List<String> aListOfUniqueSmiles);

    /**
     * Method to generate count fingerprint all by itself by comparing the unique SMILES
     *
     * @param aMapToCountMoleculeFragmentSmiles a Map containing all fragments of the molecule as unique SMILES and the
     *    corresponding frequencies of the generated fragments.
     *    Key = unique SMILES, Value = frequency
     * @return ICountFingerprint
     */
    ICountFingerprint getCountFingerprint(Map<String, Integer> aMapToCountMoleculeFragmentSmiles);

    /**
     * Method to generate count fingerprint all by itself by comparing the unique SMILES.
     *
     * @param aListToCountMoleculeFragmentSmiles a list containing all fragments of the molecule as unique SMILES.
     *    If a fragment is present more than once in the molecule, then the SMILES
     *    corresponding to the fragment is also present more than once in the list.
     * @return
     */
    ICountFingerprint getCountFingerprint(List<String> aListToCountMoleculeFragmentSmiles);



























}
