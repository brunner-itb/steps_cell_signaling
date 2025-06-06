
import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *
from steps.visual import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# mesh = TetMesh.LoadGmsh('/home/pb/steps_cell_signaling/gmsh_tutorial/' + 't1.msh')
# mesh = TetMesh.LoadAbaqus("/home/pb/steps_cell_signaling/Patrick/mesh_loading_test/Finger_0.4.inp")
mesh = TetMesh.LoadAbaqus("/home/pb/steps_cell_signaling/Patrick/meshes/elipsoid_4.5.inp")


def plotTriangles(ax, tris, color):
    ax.add_collection(Poly3DCollection(
        [tri.verts for tri in tris],
        facecolor=color,
        edgecolors='black',
        linewidth=0.1
    ))


fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(projection='3d')

plotTriangles(ax, mesh.surface, (0.5, 0.5, 0.5, 0.1))
ax.set_xlim(mesh.bbox.min.x, mesh.bbox.max.x)
ax.set_ylim(mesh.bbox.min.y, mesh.bbox.max.y)
ax.set_zlim(mesh.bbox.min.z, mesh.bbox.max.z)
ax.set_xlabel('x position [m]')
ax.set_ylabel('y position [m]')
ax.set_zlabel('z position [m]')
ax.set_title("test")
plt.show()



import os

mesh_folder = "/home/pb/steps_cell_signaling/Patrick/meshes/"
inp_files = [f for f in os.listdir(mesh_folder) if f.endswith('.inp')]

failed_ones = []

for file_name in inp_files:
    try:
        mesh = TetMesh.LoadAbaqus(mesh_folder + file_name)
        def plotTriangles(ax, tris, color):
            ax.add_collection(Poly3DCollection(
                [tri.verts for tri in tris],
                facecolor=color,
                edgecolors='black',
                linewidth=0.1
            ))

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(projection='3d')

        plotTriangles(ax, mesh.surface, (0.5, 0.5, 0.5, 0.1))
        ax.set_xlim(mesh.bbox.min.x, mesh.bbox.max.x)
        ax.set_ylim(mesh.bbox.min.y, mesh.bbox.max.y)
        ax.set_zlim(mesh.bbox.min.z, mesh.bbox.max.z)
        ax.set_xlabel('x position [m]')
        ax.set_ylabel('y position [m]')
        ax.set_zlabel('z position [m]')
        ax.set_title(file_name)
        plt.show()
    except:
        print("failed on:", file_name)
        failed_ones.append(file_name)



['AllowNestedDocks', 'AllowTabbedDocks', 'AnimatedDocks', 'DockOption', 'DockOptions', 'DrawChildren', 'DrawWindowBackground', 'ForceTabbedDocks', 'GroupedDragging', 'IgnoreMask', 'PaintDeviceMetric', 'PdmDepth', 'PdmDevicePixelRatio', 'PdmDevicePixelRatioScaled', 'PdmDpiX', 'PdmDpiY', 'PdmHeight', 'PdmHeightMM', 'PdmNumColors', 'PdmPhysicalDpiX', 'PdmPhysicalDpiY', 'PdmWidth', 'PdmWidthMM', 'RenderFlag', 'RenderFlags', 'VerticalTabs', '__class__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattr__', '__getattribute__', '__getstate__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__lt__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', '_updateView', 'acceptDrops', 'accessibleDescription', 'accessibleName', 'actionEvent', 'actions', 'activateWindow', 'addAction', 'addActions', 'addDockWidget', 'addItem', 'addItems', 'addToolBar', 'addToolBarBreak', 'adjustSize', 'autoFillBackground', 'backgroundRole', 'baseSize', 'blockSignals', 'bound_max', 'bound_min', 'center', 'centralWidget', 'changeEvent', 'childAt', 'childEvent', 'children', 'childrenRect', 'childrenRegion', 'clearFocus', 'clearMask', 'close', 'closeEvent', 'colorCount', 'connectNotify', 'contentsMargins', 'contentsRect', 'contextMenuEvent', 'contextMenuPolicy', 'corner', 'create', 'createPopupMenu', 'createWindowContainer', 'cursor', 'customContextMenuRequested', 'customEvent', 'deleteLater', 'depth', 'destroy', 'destroyed', 'devType', 'devicePixelRatio', 'devicePixelRatioF', 'devicePixelRatioFScale', 'disconnect', 'disconnectNotify', 'dockOptions', 'dockWidgetArea', 'documentMode', 'dragEnterEvent', 'dragLeaveEvent', 'dragMoveEvent', 'dropEvent', 'dumpObjectInfo', 'dumpObjectTree', 'dynamicPropertyNames', 'effectiveWinId', 'ensurePolished', 'enterEvent', 'event', 'eventFilter', 'find', 'findChild', 'findChildren', 'focusInEvent', 'focusNextChild', 'focusNextPrevChild', 'focusOutEvent', 'focusPolicy', 'focusPreviousChild', 'focusProxy', 'focusWidget', 'font', 'fontInfo', 'fontMetrics', 'foregroundRole', 'frameGeometry', 'frameSize', 'geometry', 'getContentsMargins', 'getItems', 'grab', 'grabGesture', 'grabKeyboard', 'grabMouse', 'grabShortcut', 'graphicsEffect', 'graphicsProxyWidget', 'hasFocus', 'hasHeightForWidth', 'hasMouseTracking', 'hasTabletTracking', 'height', 'heightForWidth', 'heightMM', 'hide', 'hideEvent', 'hideItem', 'iconSize', 'iconSizeChanged', 'id', 'inherits', 'initPainter', 'inputMethodEvent', 'inputMethodHints', 'inputMethodQuery', 'insertAction', 'insertActions', 'insertToolBar', 'insertToolBarBreak', 'installEventFilter', 'isActiveWindow', 'isAncestorOf', 'isAnimated', 'isDockNestingEnabled', 'isEnabled', 'isEnabledTo', 'isFullScreen', 'isHidden', 'isLeftToRight', 'isMaximized', 'isMinimized', 'isModal', 'isRightToLeft', 'isSeparator', 'isSignalConnected', 'isVisible', 'isVisibleTo', 'isWidgetType', 'isWindow', 'isWindowModified', 'isWindowType', 'keyPressEvent', 'keyReleaseEvent', 'keyboardGrabber', 'killTimer', 'layout', 'layoutDirection', 'leaveEvent', 'locale', 'logicalDpiX', 'logicalDpiY', 'lower', 'main_axis', 'mapFrom', 'mapFromGlobal', 'mapFromParent', 'mapTo', 'mapToGlobal', 'mapToParent', 'mask', 'maximumHeight', 'maximumSize', 'maximumWidth', 'menuBar', 'menuWidget', 'metaObject', 'metric', 'minimumHeight', 'minimumSize', 'minimumSizeHint', 'minimumWidth', 'mouseDoubleClickEvent', 'mouseGrabber', 'mouseMoveEvent', 'mousePressEvent', 'mouseReleaseEvent', 'move', 'moveEvent', 'moveToThread', 'nativeEvent', 'nativeParentWidget', 'nextInFocusChain', 'normalGeometry', 'num_display', 'objectName', 'objectNameChanged', 'overrideWindowFlags', 'overrideWindowState', 'paintEngine', 'paintEvent', 'paintingActive', 'palette', 'parent', 'parentWidget', 'physicalDpiX', 'physicalDpiY', 'pos', 'previousInFocusChain', 'property', 'pyqtConfigure', 'raise_', 'receivers', 'rect', 'refresh', 'releaseKeyboard', 'releaseMouse', 'releaseShortcut', 'removeAction', 'removeDockWidget', 'removeEventFilter', 'removeItem', 'removeToolBar', 'removeToolBarBreak', 'render', 'repaint', 'resize', 'resizeDocks', 'resizeEvent', 'restoreDockWidget', 'restoreGeometry', 'restoreState', 'rotate', 'rotateItem', 'rotateItems', 'saveGeometry', 'saveState', 'scale', 'screen', 'scroll', 'sender', 'senderSignalIndex', 'setAcceptDrops', 'setAccessibleDescription', 'setAccessibleName', 'setAnimated', 'setAttribute', 'setAutoFillBackground', 'setBackgroundRole', 'setBaseSize', 'setCentralWidget', 'setContentsMargins', 'setContextMenuPolicy', 'setCorner', 'setCursor', 'setDisabled', 'setDockNestingEnabled', 'setDockOptions', 'setDocumentMode', 'setEnabled', 'setFixedHeight', 'setFixedSize', 'setFixedWidth', 'setFocus', 'setFocusPolicy', 'setFocusProxy', 'setFont', 'setForegroundRole', 'setGeometry', 'setGraphicsEffect', 'setHidden', 'setIconSize', 'setInputMethodHints', 'setLayout', 'setLayoutDirection', 'setLocale', 'setMask', 'setMaximumHeight', 'setMaximumSize', 'setMaximumWidth', 'setMenuBar', 'setMenuWidget', 'setMinimumHeight', 'setMinimumSize', 'setMinimumWidth', 'setMouseTracking', 'setObjectName', 'setPalette', 'setParent', 'setProperty', 'setShortcutAutoRepeat', 'setShortcutEnabled', 'setSizeIncrement', 'setSizePolicy', 'setStatusBar', 'setStatusTip', 'setStyle', 'setStyleSheet', 'setTabOrder', 'setTabPosition', 'setTabShape', 'setTabletTracking', 'setToolButtonStyle', 'setToolTip', 'setToolTipDuration', 'setUnifiedTitleAndToolBarOnMac', 'setUpdatesEnabled', 'setVisible', 'setWhatsThis', 'setWindowFilePath', 'setWindowFlag', 'setWindowFlags', 'setWindowIcon', 'setWindowIconText', 'setWindowModality', 'setWindowModified', 'setWindowOpacity', 'setWindowRole', 'setWindowState', 'setWindowTitle', 'sharedPainter', 'show', 'showEvent', 'showFullScreen', 'showItem', 'showMaximized', 'showMinimized', 'showNormal', 'signalsBlocked', 'size', 'sizeHint', 'sizeIncrement', 'sizePolicy', 'splitDockWidget', 'stackUnder', 'startTimer', 'staticMetaObject', 'statusBar', 'statusTip', 'style', 'styleSheet', 'tabPosition', 'tabShape', 'tabifiedDockWidgetActivated', 'tabifiedDockWidgets', 'tabifyDockWidget', 'tabletEvent', 'takeCentralWidget', 'testAttribute', 'thread', 'timerEvent', 'toolBarArea', 'toolBarBreak', 'toolButtonStyle', 'toolButtonStyleChanged', 'toolTip', 'toolTipDuration', 'tr', 'underMouse', 'ungrabGesture', 'unifiedTitleAndToolBarOnMac', 'unsetCursor', 'unsetLayoutDirection', 'unsetLocale', 'update', 'updateGeometry', 'updateItems', 'updateMicroFocus', 'updatesEnabled', 'visibleRegion', 'whatsThis', 'wheelEvent', 'widget', 'width', 'widthMM', 'winId', 'window', 'windowFilePath', 'windowFlags', 'windowHandle', 'windowIcon', 'windowIconChanged', 'windowIconText', 'windowIconTextChanged', 'windowModality', 'windowOpacity', 'windowRole', 'windowState', 'windowTitle', 'windowTitleChanged', 'windowType', 'x', 'y']







