public class Boost {
    public static void main(String[] args) {
	StdAudio.play("Comet.wav");
	try {
	    Thread.sleep(2000);
	} catch (InterruptedException ex) {

	}
	System.exit(0);
    }
}