package nl.tudelft.bt.model.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;

/**
 * Implementes a class that may be used to read a simulation results file from
 * disk and get values from itr
 * 
 * @author jxavier
 */
public class SimulationResultParser {
	private ArrayList _variableList;

	private class ResultVariableSeries {
		private String _name;
		private ArrayList _values;
		/**
		 * Create a new instance
		 * 
		 * @param name
		 *            of the variable series, taken from header of
		 *            simulationsResult.txt
		 */
		private ResultVariableSeries(String name) {
			_name = name;
			_values = new ArrayList();
		}

		/**
		 * Add a new value to the list
		 * 
		 * @param v
		 */
		private void addValue(float v) {
			_values.add(new Float(v));
		}
		
		/**
		 * Retrun the float value stored last
		 * 
		 * @return the last value stored
		 */
		private float getLastValue() {
			Float f = ((Float)(_values.get(_values.size()-1)));
			return f.floatValue();
		}

		
		/* (non-Javadoc)
		 * @see java.lang.Object#equals(java.lang.Object)
		 */
		public boolean equals(Object arg0) {
			String s =  toString() + " bulk concentration [M/L^3]";
			return (arg0.toString()).equals(s);
		}
		/* (non-Javadoc)
		 * @see java.lang.Object#toString()
		 */
		public String toString() {
			return _name;
		}
	}

	/**
	 * @param fileName
	 * @param sp
	 * @param coreSpecificMass
	 * @param coreSpecies
	 *            The species constituting the core
	 * @param capsuleSpecificMass
	 * @param capsuleSpecies
	 *            The species contituting the capsule
	 */
	public SimulationResultParser(String fileName) {
		// read and parse the file
		File d = new File(fileName);
		try {
			BufferedReader br = new BufferedReader(new FileReader(d));
			String line;
			String[] tokens;
			_variableList = new ArrayList();
			// create the variable lists from the first line
			// (the simulation results header)
			line = br.readLine();
			tokens = line.split("\t");
			for (int i = 0; i < tokens.length; i++) {
				_variableList.add(new ResultVariableSeries(tokens[i]));
			}
			// read the values from the remaining lines
			while ((line = br.readLine()) != null) {
				tokens = line.split("\t");
				for (int i = 0; i < tokens.length; i++) {
					float v = Float.parseFloat(tokens[i]);
					((ResultVariableSeries) (_variableList.get(i))).addValue(v);
				}
			}
			br.close();
		} catch (FileNotFoundException e) {
			throw new ModelRuntimeException("file " + fileName + " not found.");
		} catch (IOException e) {
			throw new ModelRuntimeException(e.toString());
		}
	}

	/**
	 * Get the value for a bulk concentration from the values read from a file
	 * 
	 * @param s
	 *            dolute species to get the bulk concentration
	 * @param i
	 *            iteration at which to get it
	 * @return the value of the bulk concentration for s at iteration i
	 */
	public float getValueForBulkConcentration(SoluteSpecies s) {
		ResultVariableSeries vs = new ResultVariableSeries(s.getName());
		vs = (ResultVariableSeries) (_variableList.get(_variableList
				.indexOf(vs)));
		return vs.getLastValue();
	}
}