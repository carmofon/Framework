package nl.tudelft.bt.model.apps.output;

import java.util.ArrayList;
import java.util.Iterator;

import nl.tudelft.bt.model.Model;

/**
 * Writes the position and other properties (according to biomass particle
 * method toString) to file
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class SimulationResultsWriter extends StateWriter {
	static final public String FILENAME = "simulationResults.txt";
	static final private String DELIMITER = "\t";

	private ArrayList _variableSeries;

	/**
	 * Builds a string with the name of each of the variables that will be
	 * written to disk at each iteration
	 */
	private void initializeParametersWriting() {
		// write the table header
		String ht = new String();
		ht += Model.model().getIterationParametersHeader();
		if (_variableSeries != null)
			for (Iterator iter = _variableSeries.iterator(); iter.hasNext();) {
				VariableSeries element = (VariableSeries) iter.next();
				ht += DELIMITER + element.getYLabel();
			}
		ht += "\n";
		// write to the file
		appendToFile(FILENAME, ht);
		// Write iteration 0
		write();
	}

	/**
	 * Adds a variable series so that its value is written to file
	 * simulationParameters.txt at each iteration
	 * 
	 * @param s
	 *            variable seris to add
	 */
	public void addSeries(VariableSeries s) {
		if (_variableSeries == null)
			_variableSeries = new ArrayList();
		_variableSeries.add(s);
	}

	/*
	 * (non-Javadoc)
	 */
	public void write() {
		String vals = new String();
		vals += Model.model().getIterationParameters();
		if (_variableSeries != null)
			for (Iterator iter = _variableSeries.iterator(); iter.hasNext();) {
				VariableSeries element = (VariableSeries) iter.next();
				vals += DELIMITER + element.getLastY();
			}
		vals += "\n";
		// if the file does not exist yet, create a new one with headers
		if (!fileExists(FILENAME))
			initializeParametersWriting();
		appendToFile(FILENAME, vals);
	}
}
