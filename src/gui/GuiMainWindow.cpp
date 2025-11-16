//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-10-14 17:22:03 taubin>
//------------------------------------------------------------------------
//
// GuiMainWindow.cpp
//
// Software developed for the course
// Digital Geometry Processing
// Copyright (c) 2025, Gabriel Taubin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Brown University nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL GABRIEL TAUBIN BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <iostream>

#include "GuiStrings.hpp"
#include "GuiApplication.hpp"
#include "GuiMainWindow.hpp"
#include "GuiGLWidget.hpp"
#include "GuiAboutDialog.hpp"

// #include "GuiPanelGeneral.hpp"
// #include "GuiPanelSmoothing.hpp"

#include <QMenuBar>
#include <QGroupBox>
#include <QStatusBar>
#include <QFileDialog>
#include <QFileInfo>
#include <QRect>
#include <QMargins>

#include "io/LoaderWrl.hpp"
#include "io/SaverWrl.hpp"
#include "io/LoaderStl.hpp"
#include "io/SaverStl.hpp"
#include "io/LoaderPly.hpp"
#include "io/SaverPly.hpp"
#include "io/LoaderOff.hpp"
#include "io/SaverOff.hpp"
#include "io/LoaderObj.hpp"
#include "io/SaverObj.hpp"

int GuiMainWindow::_timerInterval = 20;

//////////////////////////////////////////////////////////////////////
// Registered Panels

GuiPanelFactoryMap GuiMainWindow::_registeredPanelFactories;
GuiPanelMap        GuiMainWindow::_registeredPanels;

GuiPanelFactoryMap& GuiMainWindow::getRegisteredPanelFactories() {
  return _registeredPanelFactories;
}
GuiPanelMap& GuiMainWindow::getRegisteredPanels() {
  return _registeredPanels;
}

bool GuiMainWindow::registerPanelFactory
(QString name, GuiPanelFactory factory) {
  // TODO Sun Oct 16 17:02:57 2016
  // _registeredPanels.find(name);
  // delete old panel ?
  getRegisteredPanelFactories().insert(name,factory);
  getRegisteredPanels().insert(name,(GuiPanel*)0);
  return true;
}

GuiPanel* GuiMainWindow::getCurrentPanel() {
  // QString name = _panelCombo->currentText();
  // std::cout << "GuiMainWindow::getCurrentPanel() {\n";
  // std::cout << "  name = \"" << qPrintable(name) << "\" \n";
  // std::cout << "}\n";
  // return (GuiPanel*)0;
  return _currentPanel;
}

void GuiMainWindow::on_panelCombo_currentIndexChanged(int index) {
  (void)index;
  setCurrentPanel(_panelCombo->currentText());
}

void GuiMainWindow::selectPanel(QString name) {
  int indx = _panelCombo->findText(name);
  if(indx>=0) {
    _panelCombo->setCurrentIndex(indx); // emits currentIndexChanged()
  }
}

void GuiMainWindow::setCurrentPanel(QString name) {
  
  // look for panel by name
  GuiPanelMap::iterator iPanel = getRegisteredPanels().find(name);
  if(iPanel == getRegisteredPanels().end()) { // name notfound
    qDebug("  panel \"%s\" not found",qPrintable(name));
  } else {
    GuiPanel* panel = iPanel.value();
    // create panel on first call
    if(panel==(GuiPanel*)0) {
      // find the corresponding factory (should not fail)
      GuiPanelFactoryMap::iterator iFactory =
        getRegisteredPanelFactories().find(name);
      GuiPanelFactory factory = iFactory.value();
      // create the panel
      panel = factory(_centralWidget);
      panel->setMinimumWidth(_toolsWidth);
    }
    if(_currentPanel!=panel) {
      // hide current panel
      if(_currentPanel!=(GuiPanel*)0)
        _currentPanel->setVisible(false);
      _currentPanel = panel;
      // show newpanel
      if(_currentPanel!=(GuiPanel*)0) {
        _currentPanel->setVisible(true);
        _currentPanel->show();
      }
    }
  }
  _refresh();
}

//////////////////////////////////////////////////////////////////////
GuiMainWindow::~GuiMainWindow() {
}

GuiMainWindow::GuiMainWindow():
  _startUpPanelName("General"),
  _currentPanel(),
  _bdryTop(5),
  _bdryBottom(5),
  _bdryLeft(5),
  _bdryCenter(5),
  _bdryRight(5),
  _toolsWidth(300)
{

    setlocale(LC_ALL, "POSIX");


  setWindowIcon(QIcon("qt.icns"));
  setWindowTitle(QString("DGP2025-A4 | Student : %1").arg(STUDENT_NAME));
  
  LoaderWrl* wrlLoader = new LoaderWrl();
  _loader.registerLoader(wrlLoader);
  SaverWrl* wrlSaver = new SaverWrl();
  _saver.registerSaver(wrlSaver);

  LoaderStl* stlLoader = new LoaderStl();
  _loader.registerLoader(stlLoader);
  SaverStl* stlSaver = new SaverStl();
  _saver.registerSaver(stlSaver);

  LoaderPly* plyLoader = new LoaderPly();
  _loader.registerLoader(plyLoader);
  SaverPly* plySaver = new SaverPly();
  _saver.registerSaver(plySaver);

  LoaderOff* offLoader = new LoaderOff();
  _loader.registerLoader(offLoader);
  SaverOff* offSaver = new SaverOff();
  _saver.registerSaver(offSaver);

  LoaderObj* objLoader = new LoaderObj();
  _loader.registerLoader(objLoader);
  SaverObj* objSaver = new SaverObj();
  _saver.registerSaver(objSaver);

  // QColor clearColor    = qRgb(175, 200, 150);
  QColor clearColor    = qRgb(200, 200, 200);
  QColor materialColor = qRgb(225, 150, 75);

  _centralWidget = new QWidget(this);
  setCentralWidget(_centralWidget);

  _glWidget    =  new GuiGLWidget(this, clearColor, materialColor);
  _glWidget->setParent(_centralWidget);

  _panelCombo = new QComboBox(_centralWidget);
  {
    QFont font;
    font.setPointSize(10);
    _panelCombo->setFont(font);
  }
  connect(_panelCombo, SIGNAL(currentIndexChanged(int)),
          this, SLOT(on_panelCombo_currentIndexChanged(int)));
      
  // fill the panel combo
  auto panelNames = getRegisteredPanels().keys();
  foreach (auto item, panelNames) {
    _panelCombo->insertItem(999, item);
  }

  _glWidget->show();

  _statusBar = new QStatusBar(this);
  {
    QFont font;
    font.setPointSize(9);
    _statusBar->setFont(font);
  }
  setStatusBar(_statusBar);
  // _statusBar->hide();

  //////////////////////////////////////////////////
  QMenu *fileMenu = menuBar()->addMenu("&File");

  QAction *exit = new QAction("E&xit" , fileMenu);
  fileMenu->addAction(exit);
  connect(exit, SIGNAL(triggered(bool)),
          this, SLOT(onMenuFileExit()));
  
  QAction *load = new QAction("Load" , fileMenu);
  fileMenu->addAction(load);
  connect(load, SIGNAL(triggered(bool)),
          this, SLOT(onMenuFileLoad()));
  
  QAction *save = new QAction("Save" , fileMenu);
  fileMenu->addAction(save);
  connect(save, SIGNAL(triggered(bool)),
          this, SLOT(onMenuFileSave()));
  
  QAction *qtLogo = new QAction("Qt Logo" , fileMenu);
  fileMenu->addAction(qtLogo);
  connect(qtLogo, SIGNAL(triggered(bool)),
          _glWidget, SLOT(setQtLogo()));

  //////////////////////////////////////////////////
  QMenu *toolsMenu = menuBar()->addMenu("&Tools");

  QAction *showTools = new QAction("Show" , toolsMenu);
  toolsMenu->addAction(showTools);
  connect(showTools, SIGNAL(triggered(bool)),
          this, SLOT(onMenuToolsShow()));

  QAction *hideTools = new QAction("Hide" , toolsMenu);
  toolsMenu->addAction(hideTools);
  connect(hideTools, SIGNAL(triggered(bool)),
          this, SLOT(onMenuToolsHide()));

  //////////////////////////////////////////////////
  QMenu *helpMenu = menuBar()->addMenu("&Help");

  QAction *about = new QAction("About QtOpenGL",helpMenu);
  helpMenu->addAction(about);
  connect(about, SIGNAL(triggered(bool)),
           this, SLOT(onMenuHelpAbout()));

  // for animation
  _timer = new QTimer(this);
  _timer->setInterval(_timerInterval);
  connect(_timer, SIGNAL(timeout()), _glWidget, SLOT(update()));

  showStatusBarMessage("");
  _glWidget->setFocus();

  setMinimumSize(500,500);
  setGLWidgetSize(600,600);
}

//////////////////////////////////////////////////////////////////////
void GuiMainWindow::showEvent(QShowEvent *event) {
  (void)event;
  selectPanel(_startUpPanelName);
}

//////////////////////////////////////////////////////////////////////
void GuiMainWindow::showStatusBarMessage(const QString & message) {
  _statusBar->showMessage(message);
}

//////////////////////////////////////////////////////////////////////
SceneGraph* GuiMainWindow::loadSceneGraph(const char* fname) {
  static char str[1024];
  snprintf(str,1024,"Trying to load \"%s\" ...",fname);
  showStatusBarMessage(QString(str));
  SceneGraph* pWrl = new SceneGraph();
  if(_loader.load(fname,*pWrl)) { // if success
    snprintf(str,1024,"Loaded \"%s\"",fname);
    pWrl->updateBBox();
    _glWidget->setSceneGraph(pWrl,true);
    // updateState();
    GuiPanel* currentPanel = getCurrentPanel();
    if(currentPanel) currentPanel->updateState();
  } else {
    snprintf(str,1024,"Unable to load \"%s\"",fname);
    delete pWrl;
    pWrl = (SceneGraph*)0;
  }
  showStatusBarMessage(QString(str));
  emit fileLoaded();
  return pWrl;
}

//////////////////////////////////////////////////////////////////////
void GuiMainWindow::onMenuFileLoad() {

  // stop animation
  _timer->stop();

  QSettings& settings = getApp()->settings();
  QString dataDir = settings.value("dataDir").toString();

  QFileDialog fileDialog(this);
  fileDialog.setFileMode(QFileDialog::ExistingFile); // allowed to select only one 
  fileDialog.setAcceptMode(QFileDialog::AcceptOpen);
  if(dataDir!="") fileDialog.setDirectory(dataDir);
  fileDialog.setNameFilter(tr("3D Files (*.wrl *.stl *.ply *.off *.obj)"));

  std::string filePath;
  if(fileDialog.exec()) {
    QStringList fileNames = fileDialog.selectedFiles();
    if(fileNames.size()>0)
      filePath = fileNames.at(0).toStdString();
  }

  QFileInfo fInfo(filePath.c_str());
  QString dir = fInfo.dir().absolutePath();
  settings.setValue("dataDir",dir);

  if (filePath.empty()) {
    showStatusBarMessage("load filePath is empty");
  } else {
    loadSceneGraph(filePath.c_str());
  } 

  // restart animation
  _timer->start(_timerInterval);
}

//////////////////////////////////////////////////////////////////////
void GuiMainWindow::onMenuFileSave() {

  // stop animation
  _timer->stop();

  // QSettings& settings = getApp()->settings();

  SceneGraph*  pWrl = _glWidget->getSceneGraph();
  std::string& url  = pWrl->getUrl(); // this is actually a file path
  QString wrlFilePath = url.c_str();
  QFileInfo fInfo(wrlFilePath);
  QString dir = fInfo.dir().absolutePath();
  QString name = fInfo.fileName();

  QFileDialog fileDialog(this);
  fileDialog.setFileMode(QFileDialog::AnyFile); // allowed to select only one 
  fileDialog.setAcceptMode(QFileDialog::AcceptSave);
  fileDialog.setNameFilter(tr("3D Files (*.wrl *.stl *.ply *.off *.obj)"));
  if(url!="") {
    // fileDialog.setDirectory(dir);
    fileDialog.selectFile(wrlFilePath);
  }

  // fileDialog.setDirectory(QString);

  std::string filePath;
  if(fileDialog.exec()) {
    QStringList fileNames = fileDialog.selectedFiles();
    if(fileNames.size()>0)
      filePath = fileNames.at(0).toStdString();
  }

  if (filePath.empty()) {
    showStatusBarMessage("save filePath is empty");
  } else {

    static char str[1024];

    snprintf(str,1024,"Saving \"%s\" ...",filePath.c_str());
    showStatusBarMessage(QString(str));

    if(_saver.save(filePath.c_str(),*pWrl)) {
      snprintf(str,1024,"Saved \"%s\"", filePath.c_str());
    } else {
      snprintf(str,1024,"Unable to save \"%s\"",filePath.c_str());
    }

    showStatusBarMessage(QString(str));
  }

  // restart animation
  _timer->start(_timerInterval);
}

//////////////////////////////////////////////////////////////////////
void GuiMainWindow::onMenuFileExit() {
  close();
}

//////////////////////////////////////////////////////////////////////
void GuiMainWindow::onMenuToolsShow() {
  int c3dWidth  = getGLWidgetWidth();
  int c3dHeight = getGLWidgetHeight();
  if(_currentPanel!=(GuiPanel*)0) _currentPanel->show();
  setGLWidgetSize(c3dWidth,c3dHeight);
}

//////////////////////////////////////////////////////////////////////
void GuiMainWindow::onMenuToolsHide() {
  int c3dWidth  = getGLWidgetWidth();
  int c3dHeight = getGLWidgetHeight();
  if(_currentPanel!=(GuiPanel*)0) _currentPanel->hide();
  setGLWidgetSize(c3dWidth,c3dHeight);
}

//////////////////////////////////////////////////////////////////////
void GuiMainWindow::onMenuHelpAbout() {
  GuiAboutDialog dialog(this);
  dialog.exec();
}

//////////////////////////////////////////////////////////////////////
int GuiMainWindow::getGLWidgetWidth() {
  return _glWidget->size().width();
}

//////////////////////////////////////////////////////////////////////
int GuiMainWindow::getGLWidgetHeight() {
  return _glWidget->size().height();
}

//////////////////////////////////////////////////////////////////////
void GuiMainWindow::_resize(int width, int height) {

  int rowH      = 30;
  int rowSpace  = 3;

  QMenuBar   *menuBar   = this->menuBar();
  QStatusBar *statusBar = this->statusBar();

  int mbH = (menuBar  )?menuBar->height()  :0;
  int sbH = (statusBar)?statusBar->height():0;

  int x0 = _bdryLeft;
  int w0 = width-_bdryLeft-_bdryRight;
  int w1 = _toolsWidth;
  int y  = _bdryTop;
  int h  = height-_bdryTop-_bdryBottom-sbH-mbH;

  if(isVisible()==false ||
     (_currentPanel!=(GuiPanel*)0 && _currentPanel->isVisible())) {
    w0 -= (w1+_bdryCenter);
    int x1 = x0+w0+_bdryCenter;
    int y1 = y;
    int h1 = h;

    // _panelCombo;
    _panelCombo->setGeometry(x1,y1,w1,rowH);
    y1 += rowH+rowSpace;
    h1 -= rowH+rowSpace;
    if(_currentPanel!=(GuiPanel*)0)
      _currentPanel->setGeometry(x1,y1,w1,h1);
  }

  _glWidget->setGeometry(x0,y,w0,h);

  // updateState();
}

//////////////////////////////////////////////////////////////////////
void GuiMainWindow::_refresh() {
  _resize(width(),height());
}

//////////////////////////////////////////////////////////////////////
void GuiMainWindow::resizeEvent(QResizeEvent* event) {
  QSize size = event->size();
  _resize(size.width(),size.height());
  if(_currentPanel) _currentPanel->updateState();
}

//////////////////////////////////////////////////////////////////////
void GuiMainWindow::setGLWidgetSize(int c3dWidth, int c3dHeight) {
  int width  = _bdryLeft + c3dWidth  +  _bdryRight;

  if(isVisible()==false ||
     (_currentPanel!=(GuiPanel*)0 && _currentPanel->isVisible()))
    width  += _bdryCenter + _toolsWidth;

  QMenuBar   *menuBar   = this->menuBar();
  QStatusBar *statusBar = this->statusBar();
  int mbH = (menuBar  )?menuBar->height()  :0;
  int sbH = (statusBar)?statusBar->height():0;
  int height = _bdryTop  + c3dHeight + _bdryBottom + sbH+mbH;
  this->resize(width,height);
}

GuiViewerData& GuiMainWindow::getData() const {
  return _glWidget->getData();
}

SceneGraph* GuiMainWindow::getSceneGraph() {
  return _glWidget->getSceneGraph();
}

void GuiMainWindow::setSceneGraph(SceneGraph* pWrl, bool resetHomeView) {
  _glWidget->setSceneGraph(pWrl,resetHomeView);
}

void GuiMainWindow::resetSceneGraph() {
  _glWidget->resetSceneGraph();
}

void GuiMainWindow::refresh() {
  _glWidget->update();
}
