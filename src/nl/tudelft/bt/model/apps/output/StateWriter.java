package nl.tudelft.bt.model.apps.output;

import java.io.File;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.ModelHandler;
import nl.tudelft.bt.model.exceptions.ModelException;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.exceptions.SystemEditViolationException;

/**
 * Interface for all classes that write the model state at each iteration
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public abstract class StateWriter implements StateWriterInterface, Serializable {
	static private File _outputDirectory;
	static private ApplicationComponent _modelHandle;

	/**
	 * Set the directory to write results to. Create new directory if direcotry
	 * does not exist or delete directory contents if directory exists.
	 * 
	 * @param d
	 *            directory name
	 * @throws IOException
	 *             if SecurityException is caught
	 */
	public static void setNewResultDirectoryOrDeleteContents(String d) throws IOException {
		File dir = new File(d);
		try {
			if (dir.isDirectory())
				deleteDirectoryContents(dir);
			dir.mkdir();
		} catch (SecurityException e) {
			throw new IOException("directory " + d
					+ " could not be created for security reasons");
		}
		_outputDirectory = dir;
	}

	/**
	 * Sets an existing directory as the output directory
	 * 
	 * @param d
	 * @throws ModelRuntimeException
	 */
	public static void reUseResultDirectory(String d)
			throws ModelRuntimeException {
		File dir = new File(d);
		if (!dir.isDirectory())
			throw new ModelRuntimeException("directory " + d
					+ " does not exist");
		_outputDirectory = dir;
	}
	/**
	 * }
	 * 
	 * @return the full name (including path) of the output directory
	 */
	public static String getOutputDirectoryName() {
		return _outputDirectory.getPath();
	}

	/**
	 * Delete the contents of a directory
	 * 
	 * @param d
	 */
	private static void deleteDirectoryContents(File d) {
		File[] files = d.listFiles();
		for (int i = 0; i < files.length; i++) {
			File f = files[i];
			if (f.isDirectory())
				deleteDirectoryContents(f);
			else
				f.delete();
		}
		d.delete();
	}

	/**
	 * Set the model Handler
	 * 
	 * @param h
	 */
	public static void setModelHandle(ModelHandler h) {
		_modelHandle = h;
	}

	/**
	 *  
	 */
	public StateWriter() throws SystemEditViolationException {
		if (Model.model() == null)
			throw new SystemEditViolationException(
					"Trying to create a System handler"
							+ " without setting model handler attribute");
		if (_outputDirectory == null)
			throw new SystemEditViolationException(
					"Trying to create a System handler"
							+ " without setting output directory");
	}

	/**
	 * Return the iteration number formated as a string to be used in file name
	 * numbering
	 * 
	 * @return formatted string of iteration number
	 */
	public String getFormattedIterationNumber() {
		return Model.model().getFormatedIterationNumber();
	}

	/**
	 * Write a set of model ooutput to a given directory
	 */
	public abstract void write() throws ModelException;

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.apps.output.StateWriterInterface#dump()
	 */
	public void dump() {
		try {
			write();
		} catch (ModelException e) {
			System.out.println(e.toString());
		}
	}

	/**
	 * @param os
	 * @throws IOException
	 */
	public static void serializeStaticState(ObjectOutputStream os)
			throws IOException {
		os.writeObject(_outputDirectory);
		os.writeObject(_modelHandle);
	}

	/**
	 * @param os
	 * @throws IOException
	 * @throws ClassNotFoundException
	 */
	public static void deserializeStaticState(ObjectInputStream os)
			throws IOException, ClassNotFoundException {
		_outputDirectory = (File) (os.readObject());
		_modelHandle = (ApplicationComponent) (os.readObject());
	}

	/**
	 * Confirm the existance of the subdirectory
	 * 
	 * @param d
	 *            the subdirectory to check
	 * @return a String with the full path fo the subdirectory including closing
	 *         "/"
	 */
	public static String confirmSubDirectoryExists(String d) {
		File dir = new File(_outputDirectory + "/" + d);
		try {
			if (!dir.isDirectory())
				dir.mkdir();
			return dir.getPath() + "/";
		} catch (SecurityException e) {
			throw new SystemEditViolationException("directory " + dir
					+ " could not be created");
		}
	}

	/**
	 * Check if a file exists in the output directory
	 * 
	 * @param file
	 *            name of file to check
	 * @return true if file exists, false otherwise
	 */
	protected static boolean fileExists(String file) {
		File f = new File(_outputDirectory + "/" + file);
		return f.exists();
	}

	/**
	 * Append text to a file
	 * 
	 * @param file
	 *            file to write to
	 * @param text
	 *            the string to write
	 */
	public static void appendToFile(String file, String text) {
		// open a filewritter to append
		try {
			java.io.FileWriter fr = new java.io.FileWriter(_outputDirectory
					+ "/" + file, true);
			fr.write(text);
			fr.close();
		} catch (IOException e) {
			throw new ModelRuntimeException("Error trying to write to file "
					+ _outputDirectory + "/" + file);
		}
	}
}