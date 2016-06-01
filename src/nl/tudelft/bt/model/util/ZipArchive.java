/* 
 * Created on Mar 14, 2004 
 * by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.util;
import java.io.*;
import java.util.zip.*;
/**
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ZipArchive {
	static final int BUFFER = 2048;
	private ZipOutputStream out;
	/**
	 *  
	 */
	public ZipArchive(String zipArchive) throws IOException {
		// create file
		FileOutputStream dest = new FileOutputStream(zipArchive);
		out = new ZipOutputStream(new BufferedOutputStream(dest));
	}
	/**
	 * Add a file to a zip archive
	 * 
	 * @param fn
	 * @param zipArchive
	 */
	public final void addToZipArchiveAndDelete(File f) throws IOException{
		BufferedInputStream origin = new BufferedInputStream(
				new FileInputStream(f), BUFFER);
		out.putNextEntry(new ZipEntry(f.getName()));
		int count;
		byte data[] = new byte[BUFFER];
		while ((count = origin.read(data, 0, BUFFER)) != -1) {
			out.write(data, 0, count);
		}
		origin.close();
		f.delete();
	}
	/**
	 * Closes the zip archive
	 * 
	 * @throws IOException
	 */
	public void close() throws IOException {
		out.close();
	}
}
