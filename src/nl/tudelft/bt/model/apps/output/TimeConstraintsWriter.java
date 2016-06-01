package nl.tudelft.bt.model.apps.output;

import java.io.File;
import java.util.ArrayList;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.timeconstraint.TimeStepConstraint;

/**
 * Writes the list of time constraits to file
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class TimeConstraintsWriter extends StateWriter {
	static final public String FILEWTHSTR = "timeConstraintsCause.txt";
	static final public String FILEWTHOUTSTR = "timeConstraints.txt";

	private File _fileWithString;
	private File _fileWithoutString;

	private ArrayList _timeConstraintSeries;

	private void initializeParametersWriting() {
		// create file with string
		appendToFile(FILEWTHSTR, TimeStepConstraint.getHeaderWithString());
		// create file without string
		appendToFile(FILEWTHOUTSTR, TimeStepConstraint.getHeaderWithoutString());
	}

	/**
	 * Adds a variable series so that its value is written to file
	 * simulationParameters.txt at each iteration
	 * 
	 * @param s
	 *            variable seris to add
	 */
	public void addSeries(VariableSeries s) {
		if (_timeConstraintSeries == null)
			_timeConstraintSeries = new ArrayList();
		_timeConstraintSeries.add(s);
	}

	/*
	 * (non-Javadoc)
	 */
	public void write() {
		// check if the files where already created
		if (!fileExists(FILEWTHSTR))
			initializeParametersWriting();
		// wite file with string
		appendToFile(FILEWTHSTR, Model.model().getLastTimeConstraint()
				.writeWithName());
		// wite file without string
		appendToFile(FILEWTHOUTSTR, Model.model().getLastTimeConstraint()
				.writeWithoutName());

	}
}