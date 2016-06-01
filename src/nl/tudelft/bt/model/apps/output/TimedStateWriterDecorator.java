package nl.tudelft.bt.model.apps.output;

import java.io.Serializable;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.ModelException;
import nl.tudelft.bt.model.exceptions.SystemEditViolationException;

/**
 * Implements a decoration of a state writer. The StateWriter decorated with
 * this class will only write its output to file when
 * Model.model().writeTimedWriters() returns true.
 */
public class TimedStateWriterDecorator
		implements
			StateWriterInterface,
			Serializable {
	StateWriter _decorated;

	/**
	 * @throws SystemEditViolationException
	 */
	public TimedStateWriterDecorator(StateWriter s)
			throws SystemEditViolationException {
		_decorated = s;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.apps.output.StateWriter#write()
	 */
	public void write() throws ModelException {
		if (Model.model().writeTimedWriters())
			_decorated.write();
	}
	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.apps.output.StateWriterInterface#dump()
	 */
	public void dump() {
		try {
			_decorated.write();
		} catch (ModelException e) {
			System.out.println(e.toString());
		}
	}
}