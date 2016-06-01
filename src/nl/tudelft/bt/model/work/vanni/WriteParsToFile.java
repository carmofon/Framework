package nl.tudelft.bt.model.work.vanni;

import java.io.*;

/**
 * @author bucciv
 *
 */
public class WriteParsToFile {
	
	String path2;
	boolean append_to_file = false;
	File tempfile;
	

	public WriteParsToFile (String file_path) {
		path2 = file_path;
		
	}
	
	public  WriteParsToFile (String file_path, boolean append_value) {
		path2 = file_path;
		append_to_file = append_value;
		//tempfile = new File(path2);
	    //tempfile.delete();
	    //System.out.println("Deleted previous file :" + path2); 
	
	}

	/**
	 * Append a string to a parameter file
	 * 
	 * @param textline the string to append
	 * @throws IOException when file is not found
	 */
	public void writeToFile(String textline) throws IOException {
		FileWriter write = new FileWriter (path2, append_to_file);
		PrintWriter print_line = new PrintWriter (write);
		
		print_line.printf("%s" + "%n", textline);
		
		print_line.close();
	}
	
	public void writeToFile(float value) throws IOException {
		Float f = new Float(value);
		writeToFile (f.toString());

	}
	
}