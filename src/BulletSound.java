public class BulletSound {
    public static void main(String[] args) {
	StdAudio.play("bullet.wav");
	StdAudio.close();
	try {
	    Thread.sleep(500);
	} catch (InterruptedException ex) {

	}
	System.exit(0);
    }
}