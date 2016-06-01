/*
 * Created on Sep 8, 2004
 */
package nl.tudelft.bt.model.util;


/**
 * Implements a semaphore class for synchronization problems
 * 
 * @author jxavier
 */
public class Semaphore {
	private boolean _closed;
	/**
	 * Create a new open semaphore
	 */
	public Semaphore() {
		_closed = false;
	}
	/**
	 * open the semaphore and notify all threads waiting to continue
	 */
	public synchronized void open() {
		_closed = false;
		notifyAll();
	}
	/**
	 * close the semaphore
	 */
	public void close() {
		_closed = true;
	}
	/**
	 * call open() if semaphore is closed, call close() otherwise
	 */
	public void switchOpenClose() {
		if (_closed)
			open();
		else
			close();
	}
	/**
	 * Make a thread wait if the semaphore is closed
	 */
	public synchronized void waitIfClosed() {
		if (_closed)
			try {
				wait();
			} catch (InterruptedException e) {
				// Do nothing
				e.printStackTrace();
			}
	}
}