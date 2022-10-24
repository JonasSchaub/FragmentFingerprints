package de.unijena.cheminf.fragment.fingerprint;

import java.util.ArrayList;
import java.util.HashMap;

public interface IFingerprint {

    ArrayList generateFragmentFingerprint(ArrayList<String> aList, HashMap<String,String> aMap);

    int getFragmentFingerprintSize();

   Object getNumberPositiveBits(ArrayList<Object> aList);

    int getIndicesPositiveBits();

    //TODO sparse vector

















}
