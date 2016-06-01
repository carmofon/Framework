package nl.tudelft.bt.model.apps.output;

import java.io.Serializable;

import nl.tudelft.bt.model.bulkconcentrations.SbrBulkConcentration;
import nl.tudelft.bt.model.exceptions.ModelException;
import nl.tudelft.bt.model.exceptions.SystemEditViolationException;

/**
 * Implements a decoration of a state writer. The StateWriter decorated with
 * this class will only write its every nCycles number o cycles of an SBR
 * (sequencing batch reactor. The SBR cycle that manages this operation must be
 * se upon creation of this decorator
 */
public class SbrFullCycleStateWriterDecorator
		implements
			StateWriterInterface,
			Serializable {
	StateWriter _decorated;
	int _nCycles;
	SbrBulkConcentration.SbrCycle _cycle;

	/**
	 * @param s
	 * @param cycle
	 * @param nCycles
	 * @throws SystemEditViolationException
	 */
	public SbrFullCycleStateWriterDecorator(StateWriter s,
			SbrBulkConcentration.SbrCycle cycle, int nCycles)
			throws SystemEditViolationException {
		//assign attributes
		_decorated = s;
		_cycle = cycle;
		_nCycles = nCycles;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.apps.output.StateWriter#write()
	 */
	public void write() throws ModelException {
		//write only if number of elapsed cycles is multple of nCycles
		if (((float) _cycle.elapsedCycles() % (float) _nCycles) == 0)
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