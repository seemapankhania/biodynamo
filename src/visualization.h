#ifndef VISUALIZATION_H_
#define VISUALIZATION_H_

#include <TEveBrowser.h>
#include <TEveGeoNode.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TSystem.h>

namespace bdm {

/**
 * Singleton class, which creates graphical user interface for biodynamo
 * simulation
 */
class ROOTViz {
 public:
  static ROOTViz &GetInstance() {
    static ROOTViz instance;  // Guaranteed to be destroyed.
    // Instantiated on first use.
    return instance;
  }

  /**
   * Creates TEveManager window, initializes members
   */
  void Init() {
    TEveManager::Create();

    geom_ = new TGeoManager("Visualization", "BioDynaMo");

    // Set number of segments for approximating circles in drawing.
    // Keep it low for better performance.
    geom_->SetNsegments(4);

    mat_empty_space_ = new TGeoMaterial("EmptySpace", 0, 0, 0);
    mat_solid_ = new TGeoMaterial("Solid", .938, 1., 10000.);
    med_empty_space_ = new TGeoMedium("Empty", 1, mat_empty_space_);
    med_solid_ = new TGeoMedium("Solid", 1, mat_solid_);

    // Another way to make top volume for world. In this way it will be
    // unbounded.
    top_ = new TGeoVolumeAssembly("WORLD");

    geom_->SetTopVolume(top_);
    geom_->SetMultiThread(true);

    // connect geom_ to eve
    TGeoNode *node = geom_->GetTopNode();
    eve_top_ = new TEveGeoTopNode(geom_, node);
    gEve->AddGlobalElement(eve_top_);
    gEve->AddElement(eve_top_);

    // number of visualized nodes inside one volume. If you exceed this number,
    // ROOT will draw nothing.
    geom_->SetMaxVisNodes(max_viz_nodes_);

    gEve->GetBrowser()->GetMainFrame()->SetWindowName(
        "Biodynamo Visualization");

    init_ = true;
  }

  /**
   * Updates GLViewer of TEveManager according to current state of ECM.
   */
  void Update(Cell<Soa> &cells) {
    if (!init_)
      throw std::runtime_error("Call GUI::getInstance().Init() first!");

    if (is_geometry_closed_)
      throw std::runtime_error(
          "Geometry is already closed! Don't call GUI::Update() after "
          "GUI::CloseGeometry()!");

    top_->ClearNodes();

    for (size_t i = 0; i < cells.size(); i++) {
      auto container = new TGeoVolumeAssembly("A");
      AddBranch(cells[i], container);
      top_->AddNode(container, top_->GetNdaughters());
    }

    gEve->FullRedraw3D();
    gSystem->ProcessEvents();

    update_ = true;
  }

  /**
   * Setter for this->MaxVizNodes
   * @param number - number of visualizable nodes per volume
   */
  void SetMaxVizNodes(int number) { this->max_viz_nodes_ = number; }

  /**
   * After building the full geometry tree, geometry must be closed.
   * Closing geometry implies checking the geometry validity, fixing shapes with
   * negative parameters (run-time shapes) building the cache manager,
   * voxelizing
   * all volumes, counting the total number of physical nodes and registering
   * the manager class to the browser.
   */
  void CloseGeometry() {
    if (!update_)
      throw std::runtime_error("Call GUI::getInstance().Update() first!");

    geom_->CloseGeometry();

    is_geometry_closed_ = true;
  }

 private:
  ROOTViz() : init_(false), update_(false), max_viz_nodes_(1e6) {}
  ROOTViz(ROOTViz const &) = delete;
  ROOTViz &operator=(ROOTViz const &) = delete;

  /**
   * Recursively adds sphere and its daughters to the container.
   */
  void AddBranch(Cell<Soa>::Self<SoaRef> cell, TGeoVolume *container) {
    AddCellToVolume(cell, container);
    // to be extended for other object
  }

  /**
   * Adds sphere to the volume, its name will be: "S%d", sphereID
   */
  void AddCellToVolume(Cell<Soa>::Self<SoaRef> cell, TGeoVolume *container) {
    char name[12];
    snprintf(name, sizeof(name), "S%zu", cell.id());

    auto radius = cell.GetDiameter() / 2;
    auto massLocation = cell.GetMassLocation();
    auto x = massLocation[0];
    auto y = massLocation[1];
    auto z = massLocation[2];
    auto position = new TGeoTranslation(x, y, z);

    auto volume = geom_->MakeSphere(name, med_solid_, 0., radius);
    volume->SetLineColor(kBlue);

    container->AddNode(volume, container->GetNdaughters(), position);
  }

  /// Geometry manager
  TGeoManager *geom_;

  /// Top volumes for TGeo and TEve (world)
  TGeoVolume *top_;
  TEveGeoTopNode *eve_top_;

  /// Eve materials and medium
  TGeoMaterial *mat_empty_space_;
  TGeoMaterial *mat_solid_;
  TGeoMedium *med_empty_space_;
  TGeoMedium *med_solid_;

  // just to ensure that methods were called in correct order
  bool init_;                // true if init_ was called
  bool update_;              // true if update was called
  bool is_geometry_closed_;  // true if geometry is already closed

  /// Max visualized shapes per volume.
  int max_viz_nodes_;
};

}  // namespace bdm

#endif  // VISUALIZATION_H_
