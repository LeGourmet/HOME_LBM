public class Bubble {
  // ----------------------------------------------------- ATTRIBUTS -----------------------------------------------------
  public float rho;
  public float volume;
  public float volumeInit;
  public int numberCells;
  public boolean deprecated;
  
  // --------------------------------------------- DESTRUCTOR / CONSTRUCTOR ----------------------------------------------  
  public Bubble(){ reset(); }
  
  // ----------------------------------------------------- FUNCTIONS -----------------------------------------------------
  public void reset(){
    this.rho = 0.f;
    this.volume = 0.f;
    this.volumeInit = 0.f;
    this.numberCells = 0;
    this.deprecated = false;
  }
}
