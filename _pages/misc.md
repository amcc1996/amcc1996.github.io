---
layout: page
permalink: /misc/
title: misc
description: Miscelaneous
nav: true
nav_order: 6
---

Here is some code

```elisp
;======================================================================
;; amcc emacs configuration
;;======================================================================
;; User information
;; ----------------
(setq user-full-name "XXXX)
(setq user-mail-address "XXXX")
;;======================================================================
;; Repositories and packages
;; -------------------------
;; Adding Melpa repository
(require 'package)
(add-to-list 'package-archives
             '("melpa" . "https://melpa.org/packages/"))
(add-to-list 'package-archives
             '("melpa-stable" . "https://stable.melpa.org/packages/"))
(when (< emacs-major-version 24)
  ;; For compatibility
  (add-to-list 'package-archives '("gnu" . "http://elpa.gnu.org/packages/")))
(package-initialize)
;; Install user defined packages
(defvar amcc/packages '(
                        ;; General autocomplete package
                        auto-complete
                        ;; Syntax checking
                        flycheck
                        ;; Git porcelain
                        magit
                        ;; Org-mode
                        org
                        org-ref
                        org-bullets
                        ox-twbs
                        ox-gfm
                        ;; M-x enhacement
                        smex
                        ;; Writing tips (passive voice, etc.)
                        writegood-mode
                        ;; Dracula theme
                        dracula-theme
                        ;; Restart emacs
                        restart-emacs
                        ;; Avy: efficient motion
                        avy
                        ;; Multiple cursos
                        multiple-cursors
                        ;; Projectile for project management features
                        projectile
                        ;; Helm
                        helm
                        ;; IDO enhancements
                        ido-vertical-mode
                        ido-completing-read+
                        ;; Company (autocompletion)
                        company
                        company-auctex
                        company-reftex
                        helm-company
                        ;; Spell checking
                        ispell
                        ;; Yasnippet: template completion
                        yasnippet
                        yasnippet-snippets
                        ;; Python features
                        elpy
                        ;; Smartparens
                        smartparens
                        ;; Search on the internet
                        engine-mode
                        ;; Bind-keys
                        bind-key
                        )
  "Default packages")
(defun amcc/packages-installed-p ()
  (cl-loop for pkg in amcc/packages
           when (not (package-installed-p pkg)) do (cl-return nil)
           finally (return t)))
;; Checking if my normal package list is installed; if not, refresh it
(unless (amcc/packages-installed-p)
  (message "%s" "Refreshing package database...")
  (package-refresh-contents)
  (dolist (pkg amcc/packages)
    (when (not (package-installed-p pkg))
      (package-install pkg))))
;;======================================================================
;; Visuals and utilities
;; ---------------------
;; Select theme (set it as safe beforehand)
(add-to-list
 'custom-theme-load-path
 "~/.emacs.d/elpa/dracula-theme-20250625.2011/")
;; Set  dracula theme as safe
(setq custom-safe-themes
      (quote
       ("2d74de1cc32d00b20b347f2d0037b945a4158004f99877630afc034a674e3ab7" default)))
(load-theme 'dracula)
;; Set default font
(add-to-list 'default-frame-alist '(font . "Inconsolata-11"))
;; Highlight current row
(global-hl-line-mode 1)
;; Proper line wrapping: do not wrap lines in middle of words
(global-visual-line-mode 1)
;; Disable fringe (used tipically to show continuation lines)
(set-fringe-mode '(0 . 0))
;; Highlight parenthesis pairs (parenthesis or expression)
(show-paren-mode 1)
(setq show-paren-style 'expression)
(setq show-paren-when-point-inside-paren t)
;; Disable initial splash screen and remove initial message in scracth
(setq inhibit-splash-screen t
      initial-scratch-message nil)
;; Show an hourglass when Emacs is busy
(setq display-hourglass t)
;; Set the number of second before showing the hourglass
(setq hourglass-delay 1)
;; Set MOUSE pointer invisble when typing
(setq make-pointer-invisible t)
;; Flashes on error
(setq visible-bell t)
(setq visible-cursor nil)
;; ;; Controls visibility of the cursor
(blink-cursor-mode 0)
;; Display battery percentage
(display-battery-mode 1)
;; Display (or not) the toolbar
(tool-bar-mode 0)
;; Show scroll bar
(scroll-bar-mode 1)
;; Show menu bar
(menu-bar-mode 1)
;; Show help when mouse is stopped over a meaningful piece of text
(tooltip-mode 1)
;; Delete region when insertion starts
(delete-selection-mode t)
;; Highlights current region
(transient-mark-mode t)
;; Enables the usage of the clipboard by commands like kill and yank
(setq x-select-enable-clipboard t)
;; Shows empty lines and trailing whitespaces
(setq-default indicate-empty-lines t)
(when (not indicate-empty-lines)
  (toggle-indicate-empty-lines))
(setq-default show-trailing-whitespace t)
;; Set tab width
;; Specifying that indentation should be done with sapces and not
;; with tabs, ensuring the file will always look the same, since
;; tabs are represented differently in different editors
(setq tab-width 4
      indent-tabs-mode nil)
;; Make indent always use spaces and not tabs
(progn
  ;; make indent commands use space only (never tab character)
  (setq-default indent-tabs-mode nil)
  ;; emacs 23.1 to 26, default to t
  ;; if indent-tabs-mode is t, it means it may use tab, resulting mixed space and tab
  )
;; Remove backup files (might be useful, but insanely annoying)
(setq make-backup-files nil)
;; Using y/n instead of yes/no to answer emacs questions
(defalias 'yes-or-no-p 'y-or-n-p)
;; Clean up the buffer
(defun untabify-buffer ()
  (interactive)
  (untabify (point-min) (point-max)))
(defun indent-buffer ()
  (interactive)
  (indent-region (point-min) (point-max)))
(defun cleanup-buffer ()
  "Perform a bunch of operations on the whitespace content of a buffer."
  (interactive)
  (indent-buffer)
  (untabify-buffer)
  (delete-trailing-whitespace))
(defun cleanup-region (beg end)
  "Remove tmux artifacts from region."
  (interactive "r")
  (dolist (re '("\\\\│\·*\n" "\W*│\·*"))
    (replace-regexp re "" nil beg end)))
;; Open recent files
(recentf-mode 1)
;; Setting the maximum size of file's list
(setq recentf-max-menu-items 25)
;; Save recent files list in intervals of 5 minutes
(run-at-time nil (* 5 60) 'recentf-save-list)
;; Start GNU Emacs's windows maximized
;; (add-to-list 'default-frame-alist '(fullscreen . maximized))
;; Shell commnad keybinding
(global-set-key (kbd "C-c x") 'shell-command)
;; Adding a new paths to load plugins
;; (add-to-list 'load-path "~/.emacs.d/lisp")
;; (add-to-list 'load-path "~/.emacs.d/lisp/helm-master")
;; (add-to-list 'load-path "~/.emacs.d/lisp/fortran-tags-master")
;; Reload buffer without confirmation
(defun revert-buffer-no-confirm ()
  (interactive)
  (revert-buffer :ignore-auto :noconfirm))
;; Set keybinding to move bewtween emacs windows (Shift + Arrow)
(windmove-default-keybindings)
;; Delay to show first tooltip
(setq tooltip-delay 0)
;; Delay to show subsquent tooltips
(setq tooltip-short-delay 0)
;; Delay to hide tooltip
(setq tooltip-hide-delay 20)
;; Don't use GTK tooltip style
(setq x-gtk-use-system-tooltips nil)
;; Set intervals between auto-saves in characters
(setq auto-save-interval 100)
;; Backup settings
(setq backup-by-copying t
      ;; Don't clobber symlinks
      backup-directory-alist '(("." . "~/.emacs.d/backups/"))
      ;; Don't litter my fs tree
      delete-old-versions t
      kept-new-versions 6
      kept-old-versions 2
      version-control t
      ;; Use versioned backups
      )
;; Directory to save the auto-saves files
(setq auto-save-file-name-transforms
      `((".*" "~/.emacs.d/backups/" t)))
;; Search for -,_, tab and newline when space in isearh-forward
(defun xah-toggle-search-whitespace ()
  "Set `search-whitespace-regexp' to nil or includes hyphen lowline tab newline.
Explanation: When in isearch (M-x `isearch-forward'), space key can also stand for other chars such as hyphen lowline tab newline. It depend on a regex. It's convenient. But sometimes you want literal. This command makes it easy to toggle.

Emacs Isearch Space Toggle
http://ergoemacs.org/emacs/emacs_isearch_space.html
Version 2019-02-22"
  (interactive)
  (if (string-equal search-whitespace-regexp nil)
      (progn
        (setq search-whitespace-regexp "[-_ \t\n]+")
        (message "Space set to hyphen lowline tab newline space"))
    (progn
      (setq search-whitespace-regexp nil)
      (message "Space set to literal."))))
;; Show line numbers
(when (version<= "26.0.50" emacs-version )
  (global-display-line-numbers-mode))
;; Show columns numbers
(column-number-mode 1)
;; Set up unicode
(prefer-coding-system       'utf-8)
(set-default-coding-systems 'utf-8)
(set-terminal-coding-system 'utf-8)
(set-keyboard-coding-system 'utf-8)
(setq default-buffer-file-coding-system 'utf-8)
(setq x-select-request-type '(UTF8_STRING COMPOUND_TEXT TEXT STRING))
;;======================================================================
;; IDO (improved navigation and completion)
;; ---
;; Starting IDO Mode; makes navigating in the file system easier
(ido-mode t)
(setq ido-enable-flex-matching t
      ido-use-virtual-buffers t)
;; Vertical IDO mode
(require 'ido-vertical-mode)
(ido-vertical-mode 1)
(setq ido-vertical-show-count t)
;; Improved IDO everywhere possible (incompatible with helm)
;; (ido-everywhere 1)
;; (require 'ido-completing-read+)
;; (ido-ubiquitous-mode 1)
;;======================================================================
;; General purpose keybindings
;; ---------------------------
;; Reload file from disk without confirmation
(global-set-key (kbd "C-o") 'revert-buffer-no-confirm)
;; Remove whitespaces and empty lines on region and buffer
(global-set-key (kbd "C-x M-t") 'cleanup-region)
(global-set-key (kbd "C-c n") 'cleanup-buffer)
;; Open recent files
(global-set-key "\C-x\ \C-r" 'recentf-open-files)
;; Keybind to newline and indent
(global-set-key (kbd "RET") 'newline-and-indent)
;; Keybinding to decrease font scale
(global-set-key (kbd "C-+") 'text-scale-increase)
(global-set-key (kbd "C--") 'text-scale-decrease)
;;======================================================================
;; smartparens
;; -----------
(require 'smartparens-config)
;; Do not use Smartparens in LaTeX
;; (setq sp-ignore-modes-list '(LaTeX-mode minibuffer-incative-mode))
(setq smartparens-global-mode t)
(add-hook 'prog-mode-hook 'turn-on-smartparens-strict-mode)
(add-hook 'markdown-mode-hook 'turn-on-smartparens-strict-mode)
;; Keybindings from https://ebzzry.io/en/emacs-pairs/
(defmacro def-pairs (pairs)
  "Define functions for pairing. PAIRS is an alist of (NAME . STRING)
conses, where NAME is the function name that will be created and
STRING is a single-character string that marks the opening character.

  (def-pairs ((paren . \"(\")
              (bracket . \"[\"))

defines the functions WRAP-WITH-PAREN and WRAP-WITH-BRACKET,
respectively."
  `(progn
     ,@(cl-loop for (key . val) in pairs
             collect
             `(defun ,(read (concat
                             "wrap-with-"
                             (prin1-to-string key)
                             "s"))
                  (&optional arg)
                (interactive "p")
                (sp-wrap-with-pair ,val)))))

(def-pairs ((paren . "(")
            (bracket . "[")
            (brace . "{")
            (single-quote . "'")
            (double-quote . "\"")
            (back-quote . "`")))

(bind-keys
 :map smartparens-mode-map
 ("C-M-a" . sp-beginning-of-sexp)
 ("C-M-e" . sp-end-of-sexp)

 ("C-<down>" . sp-down-sexp)
 ("C-<up>"   . sp-up-sexp)
 ("M-<down>" . sp-backward-down-sexp)
 ("M-<up>"   . sp-backward-up-sexp)

 ("C-M-f" . sp-forward-sexp)
 ("C-M-b" . sp-backward-sexp)

 ("C-M-n" . sp-next-sexp)
 ("C-M-p" . sp-previous-sexp)

 ("C-S-f" . sp-forward-symbol)
 ("C-S-b" . sp-backward-symbol)

 ("C-<right>" . sp-forward-slurp-sexp)
 ("M-<right>" . sp-forward-barf-sexp)
 ("C-<left>"  . sp-backward-slurp-sexp)
 ("M-<left>"  . sp-backward-barf-sexp)

 ("C-M-t" . sp-transpose-sexp)
 ("C-M-k" . sp-kill-sexp)
 ("C-k"   . sp-kill-hybrid-sexp)
 ("M-k"   . sp-backward-kill-sexp)
 ("C-M-w" . sp-copy-sexp)
 ("C-M-d" . delete-sexp)

 ("M-<backspace>" . backward-kill-word)
 ("C-<backspace>" . sp-backward-kill-word)
 ([remap sp-backward-kill-word] . backward-kill-word)

 ("M-[" . sp-backward-unwrap-sexp)
 ("M-]" . sp-unwrap-sexp)

 ("C-x C-t" . sp-transpose-hybrid-sexp)

 ("C-c ("  . wrap-with-parens)
 ("C-c ["  . wrap-with-brackets)
 ("C-c {"  . wrap-with-braces)
 ("C-c '"  . wrap-with-single-quotes)
 ("C-c \"" . wrap-with-double-quotes)
 ("C-c _"  . wrap-with-underscores)
 ("C-c `"  . wrap-with-back-quotes))
;;======================================================================
;; smex - improvement to the M-x to execute commands
;; -------------------------------------------------
(require 'smex)
(smex-initialize)
;; Redefining M-x
(global-set-key (kbd "M-x") 'smex)
;; Redefining M-x for only the major modes
(global-set-key (kbd "M-X") 'smex-major-mode-commands)
;; The old M-x
(global-set-key (kbd "C-c C-c M-x") 'execute-extended-command)
;;======================================================================
;; Avy
;; ---
(require' avy)
(global-set-key (kbd "M-s M-s") 'avy-goto-char)
(global-set-key (kbd "M-s M-d") 'avy-goto-char-2)
(global-set-key (kbd "M-g M-g") 'avy-goto-timer)
(global-set-key (kbd "M-s") 'avy-goto-word-1)
;;======================================================================
;; Multiple cursors
;; ----------------
;; NOTE: Atom keybindings clash with windmove
(require 'multiple-cursors)
(global-set-key (kbd "C-S-c C-S-c") 'mc/edit-lines)
(global-set-key (kbd "C->") 'mc/mark-next-like-this)
(global-set-key (kbd "C-<") 'mc/mark-previous-like-this)
(global-set-key (kbd "C-c C-<") 'mc/mark-all-like-this)
;;======================================================================
;; Projectile
;; ----------
(projectile-mode +1)
(define-key projectile-mode-map (kbd "s-p") 'projectile-command-map)
(define-key projectile-mode-map (kbd "C-c p") 'projectile-command-map)
;;===================================================================
;; Spell checking
;; --------------
;; Requiring the package
(require 'ispell)
;; Set dictionary path
(setenv "DICPATH" "/usr/share/hunspell")
;; Set dictionary
(setenv "DICTIONARY" "pt_PT")
;; Define the program to be used
(setq ispell-program-name "hunspell")
;; Setting all dictionaries and changes to a default value
(setq ispell-extra-args   '("-d pt_PT"))
(setq ispell-current-dictionary "pt_PT")
(setq ispell-change-dictionary "pt_PT")
(defun switch-dictionary-pt-en ()
  ;; Function to switch between dictionaries
  (interactive)
  (let* ((dict ispell-current-dictionary)
         (new (if (or (string= dict "pt_PT") (string= dict nil)) "en_GB"
                "pt_PT")))
    (ispell-change-dictionary new)
    (message "Switched dictionary from %s to %s" dict new)))
;; Keybind the function
(global-set-key (kbd "C-c d") 'switch-dictionary-pt-en)
;; Spell check buffer
(global-set-key (kbd "C-c b") 'ispell-buffer)
;; Spell check word
(global-set-key (kbd "C-c w") 'ispell-word)
;;======================================================================
;; Helm (interactive minibuffer improvement)
;; ----
(require 'helm)
;; (require 'helm-config)
;; The default "C-x c" is quite close to "C-x C-c", which quits Emacs.
;; Changed to "C-c h". Note: We must set "C-c h" globally, because we
;; cannot change `helm-command-prefix-key' once `helm-config' is loaded.
(global-set-key (kbd "C-c h") 'helm-command-prefix)
(global-unset-key (kbd "C-x c"))
;; This was not present in the suggested extended config
(global-set-key (kbd "M-x") 'helm-M-x)
;; Rebind tab to run persistent action
(define-key helm-map (kbd "<tab>") 'helm-execute-persistent-action)
;; Make TAB works in terminal
(define-key helm-map (kbd "C-i") 'helm-execute-persistent-action)
;; List actions using C-z
(define-key helm-map (kbd "C-z")  'helm-select-action)
(when (executable-find "curl")
  (setq helm-google-suggest-use-curl-p t))
;; Open helm buffer inside current window, not occupy whole other window
(setq helm-split-window-in-side-p           t
      ;; Move to end or beginning of source when reaching top or bottom of source.
      helm-move-to-line-cycle-in-source     t
      ;; Search for library in `require' and `declare-function' sexp.
      helm-ff-search-library-in-sexp        t
      ;; scroll 8 lines other window using M-<next>/M-<prior>
      helm-scroll-amount                    8
      helm-ff-file-name-history-use-recentf t)
(helm-mode 1)
;;======================================================================
;; Yasnippet
;; ---------
(require 'yasnippet)
(yas-global-mode 1)
;;======================================================================
;; Company mode
;; ------------
;; Load company and its backends to AucTeX and RefTeXs
(require 'company)
(require 'company-auctex)
(require 'company-reftex)
;; ;; Start company backend to auctex
;; (company-auctex-init)
;; Start company-mode in every mode
(add-hook 'after-init-hook 'global-company-mode)
;; Set no delay to the completion: instantaneuous completions
(setq company-idle-delay t)
;; Show suggestions after entering one character
(setq company-minimum-prefix-length 1)
;; "Circular" list
(setq company-selection-wrap-around t)
;; Use tab key to cycle through suggestions ('tng' means 'tab and go')
(company-tng-configure-default)
;; Set limit to number of items of company
(setq company-tooltip-limit 20)
;; Quick help for company mode
(add-hook 'company-mode 'company-quickhelp-mode)
(setq company-quickhelp-delay nil)
;; Making C-Ret and C-Space as selection keys in company-mode
(with-eval-after-load 'company
  (define-key company-active-map (kbd "<return>") nil)
  (define-key company-active-map (kbd "RET") nil)
  (define-key company-active-map (kbd "C-SPC") #'company-complete-selection)
  (define-key company-active-map (kbd "<C-return>") #'company-complete-selection))
;; Set compnay backends (grouped and standalone)
(setq company-backends
      (quote
       ((company-auctex-macros
         company-auctex-symbols
         company-auctex-environments)
        company-auctex-bibs
        company-auctex-labels
        ;; Default backends
        company-bbdb
        company-nxml
        company-css
        company-eclim
        company-semantic
        company-clang
        company-xcode
        company-cmake
        company-capf
        company-files
        (company-dabbrev-code
         company-gtags
         company-etags
         company-keywords
         company-files)
        company-oddmuse
        company-dabbrev)
       ))
;; Adding backend for user defined snippets
(defun company-mode/backend-with-yas (backend)
  (if (and (listp backend) (member 'company-yasnippet backend))
      backend
    (append (if (consp backend) backend (list backend))
            '(:with company-yasnippet))))
(setq company-backends (mapcar #'company-mode/backend-with-yas company-backends))
;; Narrow down completion candidates with Helm
(eval-after-load 'company
  '(progn
     (define-key company-mode-map (kbd "C-:") 'helm-company)
     (define-key company-active-map (kbd "C-:") 'helm-company)))
;;======================================================================
;; Python
;; ------
;; Enable elpy to get IDE features
(elpy-enable)
;; Enable Flycheck with python
;; (when (require 'flycheck nil t)
;;   (setq elpy-modules (delq 'elpy-module-flymake elpy-modules))
;;   (add-hook 'elpy-mode-hook 'flycheck-mode))
;; Set Python command
(setq elpy-rpc-python-command "python3")
;; Set Python interpreter
(setq python-shell-interpreter "ipython3"
      python-shell-interpreter-args "-i --simple-prompt")
;; Fallback to rgrep
(defun elpy-goto-definition-or-rgrep ()
  "Go to the definition of the symbol at point, if found. Otherwise, run `elpy-rgrep-symbol'."
  (interactive)
  (ring-insert find-tag-marker-ring (point-marker))
  (condition-case nil (elpy-goto-definition)
    (error (elpy-rgrep-symbol
            (concat "\\(def\\|class\\)\s" (thing-at-point 'symbol) "(")))))
(define-key elpy-mode-map (kbd "M-.") 'elpy-goto-definition-or-rgrep)
;; Enable full font locking of inputs in the python shell
(advice-add 'elpy-shell--insert-and-font-lock
            :around (lambda (f string face &optional no-font-lock)
                      (if (not (eq face 'comint-highlight-input))
                          (funcall f string face no-font-lock)
                        (funcall f string face t)
                        (python-shell-font-lock-post-command-hook))))
(advice-add 'comint-send-input
            :around (lambda (f &rest args)
                      (if (eq major-mode 'inferior-python-mode)
                          (cl-letf ((g (symbol-function 'add-text-properties))
                                    ((symbol-function 'add-text-properties)
                                     (lambda (start end properties &optional object)
                                       (unless (eq (nth 3 properties) 'comint-highlight-input)
                                         (funcall g start end properties object)))))
                            (apply f args))
                        (apply f args))))
;;======================================================================
;; Fortran
;; -------
;; Blink matching if's
(setq fortran-blink-matching-if t)
;; Fortran settings
(setq fortran-continuation-string "&")
(setq fortran-do-indent 2)
(setq fortran-if-indent 2)
(setq fortran-structure-indent 2)
;; Fortran 90 settings
(setq f90-do-indent 2)
(setq f90-if-indent 2)
(setq f90-type-indent 2)
(setq f90-program-indent 2)
(setq f90-continuation-indent 4)
(setq f90-smart-end 'blink)
;; Set Fortran and Fortran 90 mode for appropriate extensions
(setq auto-mode-alist
      (cons '("\\.F90$" . f90-mode) auto-mode-alist))
(setq auto-mode-alist
      (cons '("\\.pf$" . f90-mode) auto-mode-alist))
(setq auto-mode-alist
      (cons '("\\.fpp$" . f90-mode) auto-mode-alist))
(setq auto-mode-alist
      (cons '("\\.F$" . fortran-mode) auto-mode-alist))
;; Hook to flycheck
(add-hook 'f90-mode-hook 'flycheck-mode)
(add-hook 'fortran-mode-hook 'flycheck-mode)
(add-hook 'f90-mode-hook 'flycheck-pos-tip-mode)
;; Show error on tooltip
(add-hook 'fortran-mode-hook 'flycheck-pos-tip-mode)
(setq flycheck-gfortran-language-standard "f2008")
(setq flycheck-gfortran-warnings '("all" "unused"))
(setq flycheck-gfortran-args '("-Wunderflow" "-Wextra"))
;; Flycheck configuration
(add-to-list 'display-buffer-alist
             ;; Error list window configuration
             `(,(rx bos "*Flycheck errors*" eos)
               (display-buffer-reuse-window
                display-buffer-in-side-window)
               (side            . bottom)
               (reusable-frames . visible)
               (window-height   . 0.23)))
;; Open error list automatically
(add-hook 'flycheck-mode-hook 'flycheck-list-errors)
;;======================================================================
;; AUCTeX | preview-LaTeX | RefTeX
;; -------------------------------
;; Set LaTeX path
(setenv "PATH" (concat "/usr/local/texlive/2020/bin/x86_64-linux/:" (getenv "PATH")))
;; Load AucTeX
(load "auctex.el" nil t t)
;; Set default LaTeX run commnad
(setq latex-run-command "pdflatex")
;; Redefine compilation including the "shell escape" option, needed
;; use externaliaztion library from TiKz
(setq LaTeX-command-style
      '(("" "%(PDF)%(latex) -shell-escape %S%(PDFout)")))
;; Enable parse on save
(setq TeX-auto-save t)
;; Enable parse on load (in order to save information about the file)
(setq TeX-parse-self t)
;; Ask what is the master file every time I create a new .tex file,
;; since i don't specify a default one
(setq-default TeX-master nil)
;; Set that the Forward and Inverse search shall be processed by SyncTeX
(setq TeX-source-correlate-method "synctex")
(add-hook 'LaTeX-mode-hook 'TeX-source-correlate-mode)
;; Enabling Forward and Inverse search in every .tex file
(setq TeX-source-correlate-start-server t)
;; Always start Forward and Inverse search
(setq TeX-view-program-selection '((output-pdf "PDF Viewer")))
;; Tell AucTeX what program to use under each circumstance
(setq TeX-view-program-list
      ;; Define viewing command for each viwer
      '(("PDF Viewer" "okular --unique %o#src:%n%b")))
(when (require 'latex nil t)
  ;; Make Forward Search work with okular (from
  ;; https://tex.stackexchange.com/questions/116633/forward-search-in-emacs-for-okular-not-working)
  (push '("%(masterdir)" (lambda nil
                           (file-truename (TeX-master-directory))))
        TeX-expand-list)
  (push '("Okular" "okular --unique %o#src:%n%(masterdir)./%b")
        TeX-view-program-list)
  (push '(output-pdf "Okular") TeX-view-program-selection))
'(TeX-view-program-list (quote (("Okular" ("okular %o") "okular --unique %o#src:%n%b"))))
;; Use PDFLaTeX for all files
(with-eval-after-load "tex" (TeX-global-PDF-mode 1))
;; Keybinding to close enviroments
(global-set-key (kbd "C-c e g") 'LaTeX-close-environment)
;; Every extra option is safe
(put 'TeX-command-extra-options 'safe-local-variable (lambda (xx) t))
;; Relativa path in images, to use with helm
(setq LaTeX-includegraphics-read-file
      'LaTeX-includegraphics-read-file-relative)
;;
;; Folding
;; -------
;; ;; I am not very fond of AucTeX folding mode
;; ;; It makes writing clumsy and not bovious at all
;; (add-hook 'LaTeX-mode-hook
;;           (lambda ()
;;             ;; Turn fold-mode with LaTeX
;;             (TeX-fold-mode 1)
;;             ;; Fold .tex files on load
;;             (add-hook 'find-file-hook 'TeX-fold-buffer t t)
;;             (add-hook 'after-change-functions
;;                       ;; Update buffer folding every time } or $ is inserted
;;                       (lambda (start end oldlen)
;;                         (when (= (- end start) 1)
;;                           (let ((char-point
;;                                  (buffer-substring-no-properties
;;                                   start end)))
;;                             (when (or (string= char-point "}")
;;                                       (string= char-point "$"))
;;                               (TeX-fold-paragraph)))))
;;                       t t)))
;; ;; Fold macros inserted with C-c C-m
;; (setq TeX-fold-auto t)
;; ;; Don't fold comments
;; (setq TeX-fold-preserve-comments t)
;;
;; Math Mode
;; ---------
;; Set new key to insert math symbols
(setq LaTeX-math-abbrev-prefix (kbd "«"))
;; Start math mode alongside with TeX
(add-hook 'LaTeX-mode-hook 'LaTeX-math-mode)
;; ;; Automatically insert braces after _ and ^
(setq TeX-electric-sub-and-superscript t)
;; Automaticallly inser \right* after a \left*
(setq LaTeX-electric-left-right-brace t)
(setq TeX-electric-math t)
(add-hook 'plain-TeX-mode-hook
          (lambda () (set (make-variable-buffer-local 'TeX-electric-math)
                          (cons "$" "$"))))
(add-hook 'LaTeX-mode-hook
          (lambda () (set (make-variable-buffer-local 'TeX-electric-math)
                          (cons "$" "$"))))
;; Add new math environments
(add-hook 'LaTeX-mode-hook 'add-my-latex-environments)
(defun add-my-latex-environments ()
  (LaTeX-add-environments
   '("equation*" LaTeX-env-label)
   '("align" LaTeX-env-label)
   )
  )
(setq font-latex-math-environments (quote
     ("equation" "eqnarray" "align" "alignat")))
;;
;; Outline Mode
;; ------------
(add-hook 'LaTeX-mode-hook 'outline-minor-mode)
;; Start outline-minor-mode in LaTeX documents
(global-set-key (kbd "C-c C-o C-h") 'hide-entry)
(global-set-key (kbd "C-c C-o C-t") 'hide-body)
(global-set-key (kbd "C-c C-o C-a") 'show-all)
(global-set-key (kbd "C-c C-o C-s") 'show-entry)
(global-set-key (kbd "C-c C-o C-n") 'outline-next-visible-heading)
(global-set-key (kbd "C-c C-o C-p") 'outline-previous-visible-heading)
(global-set-key (kbd "C-c C-o C-k") 'show-children)
;;
;; RefTeX
;; ------
(require 'reftex)
;; Adds functionality to AUCTeX labels and so on
(setq reftex-plug-into-AUCTeX t)
;; ;; Make RefTeX faster for large douments http://www.mssl.ucl.ac.uk/swift/om/sw/help/man/reftex.html
(setq reftex-enable-partial-scans t)
(setq reftex-save-parse-info t)
(setq reftex-use-multiple-selection-buffers t)
;; Start RefTeX with AUCTeX mode
(add-hook 'LaTeX-mode-hook 'turn-on-reftex)
;; ;; Start RefTeX with Emacs latex mode
(add-hook 'latex-mode-hook 'turn-on-reftex)
;; Some handy keybindings - C-c r for RefTeX and the first letter of the macro
(global-set-key (kbd "C-c r l") 'reftex-label)
(global-set-key (kbd "C-c r t") 'reftex-toc)
(global-set-key (kbd "C-c r r") 'reftex-reference)
(global-set-key (kbd "C-c r c") 'reftex-citation)
;; Smart/automatic labels
(setq  reftex-insert-label-flags '("s" t))
;; Prompt for empty optional arguments in cite macros.
(setq reftex-cite-prompt-optional-args t)
;; So that RefTeX in Org-mode knows bibliography
(setq reftex-bibliography-commands '("bibliography"
                                     "nobibliography"
                                     "addbibresource"))
;; Recognize \subcaptions, e.g. reftex-citatio
(setcdr (assoc 'caption reftex-default-context-regexps)
        "\\\\\\(rot\\|sub\\)?caption\\*?[[{]")
;; Make cref work with reftex
(setq reftex-label-alist '(AMSTeX))
(eval-after-load
    "latex"
  '(TeX-add-style-hook
    "cleveref"
    (lambda ()
      (if (boundp 'reftex-ref-style-alist)
      (add-to-list
       'reftex-ref-style-alist
       '("Cleveref" "cleveref"
         (("\\cref" ?c) ("\\Cref" ?C) ("\\cpageref" ?d) ("\\Cpageref" ?D)))))
      (reftex-ref-style-activate "Cleveref")
      (TeX-add-symbols
       '("cref" TeX-arg-ref)
       '("Cref" TeX-arg-ref)
       '("cpageref" TeX-arg-ref)
       '("Cpageref" TeX-arg-ref)))))
(setq LaTeX-reftex-cite-format-auto-activate nil)
(setq reftex-cite-format
      '((?\C-m . "\\cite[][]{%l}")
        (?s    . "\\cites[][]{%l}")
        (?p    . "\\parencites[][]{%l}")
        (?P    . "\\Parencites[][]{%l}")
        (?t . "\\textcite[]{%l}")
        (?a . "\\autocite[]{%l}")
        (?c . "\\cite[]{%l}")
        (?s . "\\smartcite[]{%l}")
        (?f . "\\footcite[]{%l}")
        (?n . "\\nocite{%l}")
        (?b . "\\blockcquote[]{%l}{}")))
;; Default bibliography file
(setq reftex-default-bibliography '("/home/amcc/Bibliography/References.bib"))
(setq bibtex-completion-bibliography "/home/amcc/Bibliography/References.bib")
;;
;; Flyspell
;; --------
;; Turn on Flyspell with LaTeX
(add-hook 'LaTeX-mode-hook (lambda () (flyspell-mode t) ))
                                        ; Spell Checking removal from some environments
(put 'LaTeX-mode 'flyspell-mode-predicate 'auctex-mode-flyspell-skip-myenv)
(defun auctex-mode-flyspell-skip-myenv ()
  (save-excursion
    (widen)
    (let ((p (point))
          (count 0))
      (not (or (and (re-search-backward "\\\\begin{\\(tikzpicture\\|circuitikz\\|forest\\)}" nil t)
                    (> p (point))
                    (or (not (re-search-forward "^\\\\end{\\(tikzpicture\\|circuitikz\\|forest\\)}" nil t))
                        (< p (point))))
               (eq 1 (progn (while (re-search-backward "`" (line-beginning-position) t)
                              (setq count (1+ count)))
                            (- count (* 2 (/ count 2)))))))))
  )
(add-hook 'LaTeX-mode-hook (lambda () (setq flyspell-generic-check-word-predicate
                                            'auctex-mode-flyspell-skip-myenv)))
;;
;; Preview LaTeX Environments
;; --------------------------
(eval-after-load "preview"
  '(add-to-list 'preview-default-preamble "\\PreviewEnvironment{tikzpicture}" t)
  )
(eval-after-load "preview" ;; might not have the expected behaviour
  '(add-to-list 'preview-default-preamble "\\PreviewEnvironment{axis}" t)
  )
(eval-after-load "preview"
  '(add-to-list 'preview-default-preamble "\\PreviewEnvironment{groupplot}" t)
  )
(eval-after-load "preview"
  '(add-to-list 'preview-default-preamble "\\PreviewEnvironment{circuitikz}" t)
  )
(eval-after-load "preview"
  '(add-to-list 'preview-default-preamble "\\PreviewEnvironment{tabular}" t)
  )
(eval-after-load "preview"
  '(add-to-list 'preview-default-preamble "\\PreviewEnvironment{figure}" t)
  )
(eval-after-load "preview"
  '(add-to-list 'preview-default-preamble "\\PreviewEnvironment{table}" t)
  )
;;
;; Fontification
;; http://lists.gnu.org/archive/html/emacs-orgmode/2009-05/msg00236.html
;; http://www.gnu.org/software/auctex/manual/auctex/Fontification-of-macros.html
(setq font-latex-match-reference-keywords
      '(
        ;; biblatex
        ("printbibliography" "[{")
        ("addbibresource" "[{")
        ;; Standard commands
        ;; ("cite" "[{")
        ("Cite" "[{")
        ("parencite" "[{")
        ("Parencite" "[{")
        ("footcite" "[{")
        ("footcitetext" "[{")
        ;; Style-specific commands
        ("textcite" "[{")
        ("Textcite" "[{")
        ("smartcite" "[{")
        ("Smartcite" "[{")
        ("cite*" "[{")
        ("parencite*" "[{")
        ("supercite" "[{")
        ;; Qualified citation lists
        ("cites" "[{")
        ("Cites" "[{")
        ("parencites" "[{")
        ("Parencites" "[{")
        ("footcites" "[{")
        ("footcitetexts" "[{")
        ("smartcites" "[{")
        ("Smartcites" "[{")
        ("textcites" "[{")
        ("Textcites" "[{")
        ("supercites" "[{")
        ;; Style-independent commands
        ("autocite" "[{")
        ("Autocite" "[{")
        ("autocite*" "[{")
        ("Autocite*" "[{")
        ("autocites" "[{")
        ("Autocites" "[{")
        ;; My custom cite commands
        ("posscite" "[{")
        ("Posscite" "[{")
        ("posscites" "[{")
        ("Posscites" "[{")
        ;; Text commands
        ("citeauthor" "[{")
        ("Citeauthor" "[{")
        ("citetitle" "[{")
        ("citetitle*" "[{")
        ("citeyear" "[{")
        ("citedate" "[{")
        ("citeurl" "[{")
        ;; Special commands
        ("fullcite" "[{")
        ;; cleveref
        ("cref" "{")
        ("Cref" "{")
        ("cpageref" "{")
        ("Cpageref" "{")
        ("cpagerefrange" "{")
        ("Cpagerefrange" "{")
        ("crefrange" "{")
        ("Crefrange" "{")
        ("labelcref" "{")))

(setq font-latex-match-textual-keywords
      '(
        ;; biblatex brackets
        ("parentext" "{")
        ("brackettext" "{")
        ("hybridblockquote" "[{")
        ;; Auxiliary Commands
        ("textelp" "{")
        ("textelp*" "{")
        ("textins" "{")
        ("textins*" "{")
        ;; subcaption
        ("subcaption" "[{")))

(setq font-latex-match-variable-keywords
      '(
        ;; amsmath
        ("numberwithin" "{")
        ;; enumitem
        ("setlist" "[{")
        ("setlist*" "[{")
        ("newlist" "{")
        ("renewlist" "{")
        ("setlistdepth" "{")
        ("restartlist" "{")
        ("crefname" "{")))
;;======================================================================
;; Org-mode
;; --------
(require 'ox-gfm)
(require 'org-bullets)
(require 'org-ref)
(require 'ox-latex)
(require 'ox-beamer)
(require 'ox-twbs)
;; Pretty bullets
(add-hook 'org-mode-hook (lambda () (org-bullets-mode 1)))
;; LaTeX rendering in situ
(setq org-latex-create-formula-image-program 'dvipng)
;; LaTeX rendering size
(setq org-format-latex-options (plist-put org-format-latex-options :scale 1.4))
;; Disable company in org-mode
(setq company-global-modes '(not org-mode))
;; Add Python interpreter to org-mode
(org-babel-do-load-languages
 'org-babel-load-languages
 '((python . t)))
(setq org-babel-python-command "python3")
;; Change hidden content symbol
(setq org-ellipsis "⤵")
;; Fontify code in code blocks
(setq org-src-fontify-natively t)
(setq org-src-tab-acts-natively t)
;; Org-mode agenda
(define-key global-map "\C-ca" 'org-agenda)
(define-key global-map "\C-ct" 'org-todo-list)
(define-key global-map "\C-cc" 'org-capture)
(define-key global-map "\C-cg" '(lambda () (interactive)(find-file "~/org/mygtd.org")))
;; Agenda customisation
(setq org-agenda-repeating-timestamp-show-all t)
(setq org-agenda-skip-deadline-if-done t)
(setq org-agenda-skip-scheduled-if-done t)
(setq org-agenda-start-on-weekday nil)
(setq org-deadline-warning-days 14)
(setq org-agenda-ndays 7)
(setq org-agenda-show-all-dates t)
(setq org-agenda-custom-commands
      '(
        ("H" "Work list"
         ((agenda)
          (tags-todo "URGENT")
          (tags-todo "PERSONAL")
          (tags-todo "READING")
          (tags-todo "WRITING")
          (tags-todo "PROGRAMMING")
          (tags-todo "TEACHING")
          (tags-todo "LINKS")))
        ("D" "Daily Action List"
         (
          (agenda "" ((org-agenda-ndays 1)
                      (org-agenda-sorting-strategy
                       (quote ((agenda time-up priority-down tag-up) )))
                      (org-deadline-warning-days 0)
                      ))))
        )
      )
;; Show closing time
(setq org-log-done t)
(setq org-closed-keep-when-no-todo t)
;; Default bibliography
(setq org-ref-default-bibliography '("/home/amcc/Bibliography/time-References.bib"))
;; TODO sequence
(setq org-todo-keywords
      '((sequence "TODO(t)" "STARTED(s)" "NEXT(n)" "WAITING(w)" "|" "DONE(d)" "CANCELLED(c)" "DEFERRED(f)")))

(setq org-todo-keyword-faces
      '(("STARTED" . (:foreground "turquoise" :weight bold))
        ("NEXT" . (:foreground "yellow" :weight bold))
        ("WAITING" . (:foreground "HotPink" :weight bold))
        ("CANCELLED" . (:foreground "red1"))
        ("DEFERRED" . (:foreground "PaleGreen"))
        ))
;; TAGS sequence
(setq org-tag-alist
      '(
        ("URGENT" . ?u)
        ("TEACHING" . ?t)
        ("READING" . ?r)
        ("WRITING" . ?w)
        ("PROGRAMMING" . ?p)
        ("LINKS" . ?l)
        ("HOME" . ?h)
        ("PERSONAL" . ?m)
        ("LEGATTO" . ?c)
        ("JIUJITSU" . ?j)
        )
      )
(setq org-tag-faces
      '(
        ("URGENT" . (:foreground "magenta" :weight bold))
        ("TEACHING" . (:foreground "orchid" :weight bold))
        ("READING" . (:foreground "MediumSeaGreen" :weight bold))
        ("WRITING" . (:foreground "SteelBlue" :weight bold))
        ("PROGRAMMING" . (:foreground "DeepPink" :weight bold))
        ("LINKS" . (:foreground "DodgerBlue" :weight bold))
        ("HOME" . (:foreground "Marron" :weight bold))
        ("PERSONAL" . (:foreground "chocolate" :weight bold))
        ("LEGATTO" . (:foreground "LightGoldenrod" :weight bold))
        ("JIUJITSU" . (:foreground "purple" :weight bold))
        )
      )
;; Directories
(setq org-directory (expand-file-name "~/org"))
(setq org-default-notes-file (concat org-directory "/mygtd.org"))
(setq org-agenda-files '("~/org"))
;; Wrap lines on start
(setq org-startup-truncated nil)
;; Skip regions for spell checking
(add-to-list 'ispell-skip-region-alist '(":\\(PROPERTIES\\|LOGBOOK\\):" . ":END:"))
(add-to-list 'ispell-skip-region-alist '("#\\+BEGIN_SRC" . "#\\+END_SRC"))
(add-to-list 'ispell-skip-region-alist '("#\\+BEGIN_EXAMPLE" . "#\\+END_EXAMPLE"))
;; Fast selection
(setq org-fast-tag-selection-single-key t)
(setq org-use-fast-todo-selection t)
;; Org-capture
(setq org-capture-templates
      '(
        ("w" "Work" entry (file+headline "~/org/mygtd.org" "Work")
         "** TODO %?\nAdded: %U\n" :prepend nil :kill-buffer t)
        ("n" "Next task" entry (file+headline "~/org/mygtd.org" "Work")
         "** NEXT %?\nAdded: %U\n" :prepend nil :kill-buffer t)
        ("p" "Personal" entry (file+headline "~/org/mygtd.org" "Personal")
         "** TODO %?\nAdded: %U\n" :prepend nil :kill-buffer t)
        ("i" "Idea" entry (file+headline "~/org/someday.org" "Ideas")
         "** %?\nAdded" :prepend nil :kill-buffer t)
        ("m" "Meeting" entry (file+headline "~/org/mygtd.org" "Calendar")
         "** TODO %? \n SCHEDULED: %^T \n" :prepend nil :kill-buffer t)
        ("j" "Note" entry (file+headline "~/org/journal.org" "Notes")
         "** %?\nAdded: %U\n" :prepend nil :kill-buffer t)
        )
      )
;; Properly export quotes and em-dashes
(setq org-export-with-smart-quotes t)
;; Configure refiling for Getting Things Done
(setq org-refile-targets (quote (
                                 ("mygtd.org" :maxlevel . 1)
                                 ("someday.org" :level . 2)
                                 )))
;; Set default column view headings: Task Total-Time Time-Stamp
(setq org-columns-default-format "%50ITEM(Task) %10CLOCKSUM %16TIMESTAMP_IA")
;;======================================================================
;; Engine mode
;; -----------
(require 'engine-mode)
(engine-mode t)
;; Set browser
(setq engine/browser-function 'browse-url-firefox)
;; Set keybinding
(engine/set-keymap-prefix (kbd "C-c s"))
;; Engines
(defengine amazon
  "http://www.amazon.com/s/ref=nb_sb_noss?url=search-alias%3Daps&field-keywords=%s")

(defengine duckduckgo
  "https://duckduckgo.com/?q=%s"
  :keybinding "d")

(defengine github
  "https://github.com/search?ref=simplesearch&q=%s")

(defengine google
  "http://www.google.com/search?ie=utf-8&oe=utf-8&q=%s"
  :keybinding "g")

(defengine google-images
  "http://www.google.com/images?hl=en&source=hp&biw=1440&bih=795&gbv=2&aq=f&aqi=&aql=&oq=&q=%s")

(defengine google-maps
  "http://maps.google.com/maps?q=%s"
  :docstring "Mappin' it up.")

(defengine project-gutenberg
  "http://www.gutenberg.org/ebooks/search/?query=%s")

(defengine qwant
  "https://www.qwant.com/?q=%s")

(defengine rfcs
  "http://pretty-rfc.herokuapp.com/search?q=%s")

(defengine stack-overflow
  "https://stackoverflow.com/search?q=%s")

(defengine twitter
  "https://twitter.com/search?q=%s")

(defengine wikipedia
  "http://www.wikipedia.org/search-redirect.php?language=en&go=Go&search=%s"
  :keybinding "w"
  :docstring "Searchin' the wikis.")

(defengine wiktionary
  "https://www.wikipedia.org/search-redirect.php?family=wiktionary&language=en&go=Go&search=%s")

(defengine wolfram-alpha
  "http://www.wolframalpha.com/input/?i=%s")

(defengine youtube
  "http://www.youtube.com/results?aq=f&oq=&search_query=%s")

(defengine powerthesaurus
  "https://www.powerthesaurus.org/%s/synonyms"
  :keybinding "p")
```