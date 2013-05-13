public class MissileSound {
    public static void main(String[] args) {
	StdAudio.play("Missile.wav");
	StdAudio.close();
	try {
	    Thread.sleep(4000);
	} catch (InterruptedException ex) {

	}
	System.exit(0);
    }
}